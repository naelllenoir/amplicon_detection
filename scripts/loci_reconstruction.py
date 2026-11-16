#!/usr/bin/env python3
"""
06/07/2026 | 09/07/2026 - Naell LENOIR
nhmmer_reconstruct_loci_featurecoord.py
====================

Usage
-----
    python3 nhmmer_reconstruct_loci.py input.tblout -o loci.gff3 [--coord-source ali]

Expected input columns from input.tblout (format nhmmer --tblout):
    target_name target_acc query_name query_acc hmmfrom hmmto \
    alifrom alito envfrom envto sqlen strand evalue score bias

Hits should have been previously ordered as follows:
    grep -v '#' table_nhmmer-test_TSPY2-Hsap_Y-Hsap.txt | sed -e 's/[[:space:]]\\+/,/g' -e 's/,/\t/g' | awk '{ if ($12=="-") { tmp=$7; $7=$8; $8=tmp; tmp2=$9; $9=$10; $10=tmp2 } print }' OFS='\t' | sort -k9,9n | cut -f1-15 >> output_orderedCoordStrand_nodesc.txt

"""

import argparse
import sys
from dataclasses import dataclass, field
from typing import List, Optional

@dataclass
class Hit:
    """Defines one nhmmer hit when parsed."""

    target: str
    strand: str
    hmmfrom: int
    hmmto: int
    alifrom: int
    alito: int
    envfrom: int
    envto: int
    sqlen: int
    evalue: float
    score: float
    bias: float

    # Coordinates to be defined
    gstart: int = field(init=False)
    gend: int = field(init=False)

    source_lines: List[int] = field(default_factory=list)

    # For coordinates, min and max are verified between from and to in case the input file has not been corrected for coordinates given the strand of the hit
    def __post_init__(self):
        self.gstart = min(self.alifrom, self.alito)
        self.gend = max(self.alifrom, self.alito)

    def set_genomic_coords(self, coord_source: str):
        if coord_source == "ali":
            # env = envelope coordinates (wider).
            self.gstart = min(self.alifrom, self.alito)
            self.gend = max(self.alifrom, self.alito)
        elif coord_source == "env":
            # ali = alignment core coordinates (narrower).
            self.gstart = min(self.envfrom, self.envto)
            self.gend = max(self.envfrom, self.envto)
        else:
            raise ValueError(f"coord_source unknown: {coord_source}")


def parse_nhmmer_tblout(path: str, coord_source: str = "ali") -> List[Hit]:
    """
    Parse nhmmer --tblout file after being ordered by env coordinates and description column removed, then append to hits list.
    """
    hits = []
    with open(path) as hit_file:
        for linenb, raw in enumerate(hit_file, start=1):
            line = raw.strip()
            if not line or line.startswith("#"):
                 continue
            fields = line.split()
            try:
                (
                    target, target_acc, query, query_acc,
                    hmmfrom, hmmto, alifrom, alito, envfrom, envto,
                    sqlen, strand, evalue, score, bias,
                ) = fields[:15]

                hit = Hit(
                    target=target,
                    strand=strand,
                    hmmfrom=int(hmmfrom),
                    hmmto=int(hmmto),
                    alifrom=int(alifrom),
                    alito=int(alito),
                    envfrom=int(envfrom),
                    envto=int(envto),
                    sqlen=int(sqlen),
                    evalue=float(evalue),
                    score=float(score),
                    bias=float(bias),
                )
            except (ValueError, TypeError) as exc:
                print(
                    f"[warn] line {linenb} ignored (parsing error: {exc})",
                    file=sys.stderr,
                )
                continue

            hit.set_genomic_coords(coord_source)
            hit.source_lines = [linenb]
            hits.append(hit)
    return hits


def merge_overlapping_hits(hits: List[Hit]) -> List[Hit]:
    """
    Fuse overlapping hits and keep wider coordinates (both genomic and hmm) with best evalue and score.
    Although it won't change locus reconstruction based on is_continuation, it will keep a correct estimation of the number of 'exon' per locus.
    """
    by_strand = {}
    for hit in hits:
        by_strand.setdefault((hit.target, hit.strand), []).append(hit)

    merged_all: List[Hit] = []
    for (target, strand), strandhits in by_strand.items():
        strandhits.sort(key=lambda hit: hit.gstart)

        merged_strandhits: List[Hit] = []
        for h in strandhits:
            if merged_strandhits and h.gstart <= merged_strandhits[-1].gend:
                prev = merged_strandhits[-1]
                prev.gend = max(prev.gend, h.gend)
                prev.gstart = min(prev.gstart, h.gstart)
                prev.hmmfrom = min(prev.hmmfrom, h.hmmfrom)
                prev.hmmto = max(prev.hmmto, h.hmmto)
                prev.evalue = min(prev.evalue, h.evalue)
                prev.score = max(prev.score, h.score)
                prev.source_lines.extend(h.source_lines)
            else:
                merged_strandhits.append(h)
        merged_all.extend(merged_strandhits)
    return merged_all


def is_continuation(prev: Hit, current: Hit, strand: str, max_size: Optional[int]) -> bool:
    """
    Test whether current hit can be the continuation of previous hit.
    """
    if max_size is not None and (current.gstart - prev.gend) > max_size:
        return False
    if strand == "+":
        return current.hmmfrom > prev.hmmfrom and current.hmmto > prev.hmmto
    elif strand == "-":
        return current.hmmfrom < prev.hmmfrom and current.hmmto < prev.hmmto
    else:
        return False


@dataclass
class Locus:
    """ Defines a reconstructed locus, e.g. a group of exons (Hit)."""

    target: str
    strand: str
    hits: List[Hit]

    @property
    def sorted_hits(self) -> List[Hit]:
        return sorted(self.hits, key=lambda hit: hit.gstart)

    @property
    def transcript_ordered_hits(self) -> List[Hit]:
        return sorted(self.hits, key=lambda hit: hit.hmmfrom)

    @property
    def start(self) -> int:
        return min(hit.gstart for hit in self.hits)

    @property
    def end(self) -> int:
        return max(hit.gend for hit in self.hits)

    @property
    def size(self) -> int:
        return self.end - self.start + 1

    @property
    def n_hits(self) -> int:
        return len(self.hits)

    @property
    def evalues(self) -> List[float]:
        return [hit.evalue for hit in self.transcript_ordered_hits]

    @property
    def scores(self) -> List[float]:
        return [hit.score for hit in self.transcript_ordered_hits]

    def format_hits(self) -> str:
        parts = [
            f"{hit.gstart}-{hit.gend}({hit.hmmfrom}-{hit.hmmto})"
            for hit in self.transcript_ordered_hits
        ]
        return ",".join(parts)

    def span_start(self, span: Optional[int] = None) -> str:
        if span is not None and (self.size) < span:
            startspan = span - self.size
        else:
            startspan = 1000
        coord_start = self.start - startspan
        return coord_start

    def span_end(self, span: Optional[int] = None) -> str:
        if span is not None and (self.size) < span:
            endspan = span - self.size
        else:
            endspan = 1000
        coord_end = self.end + endspan
        return coord_end


def build_loci(hits: List[Hit], max_size: Optional[int] = None) -> List[Locus]:
    """
    Sort hits by strand and genomic coordinates then build loci with is_continuation.
    """
    by_strand = {}
    for hit in hits:
        by_strand.setdefault((hit.target, hit.strand), []).append(hit)

    loci: List[Locus] = []
    for (target, strand), strandhits in by_strand.items():
        strandhits.sort(key=lambda hit: hit.gstart)
        current: List[Hit] = []
        for hit in strandhits:
            if current and is_continuation(current[-1], hit, strand, max_size):
                current.append(hit)
            else:
                if current:
                    loci.append(Locus(target=target, strand=strand, hits=current))
                current = [hit]
        if current:
            loci.append(Locus(target=target, strand=strand, hits=current))

    return loci



def conditional_expansion(loci: List[Locus], expansion: Optional[int] = None) -> List[Locus]:

    if expansion is None:
        return loci

    sorted_loci = sorted(loci, key=lambda locus: (locus.target, locus.strand, locus.start))
    for locus in sorted_loci:
        locus.extstart = locus.start - expansion
        locus.extend = locus.end + expansion

    for prev, current in zip(sorted_loci, sorted_loci[1:]):
        if prev.target != current.target:
            continue
        if prev.strand != current.strand:
            continue
        if prev.extend < current.extstart:
            continue
        elif prev.extend > current.extstart:
            split = (prev.end + current.start) // 2
            prev.extend = split
            current.extstart = split + 1

    return sorted_loci



def write_loci_gff3(loci: List[Locus], path: str, source: str = "nhmmer_reconstruction", span: Optional[int] = None, expansion: Optional[int] = None):
    """
    Write reconstructed loci in gff3 file.
    """
    with open(path, "w") as gff3:
        gff3.write("##gff-version 3\n")

        for i, locus in enumerate(
            sorted(loci, key=lambda locus: (locus.target, locus.start)), start=1):

            locus_id = f"locus_{i}"
            evalues_str = ",".join(f"{e:.3g}" for e in locus.evalues)
            scores_str = ",".join(f"{s:.1f}" for s in locus.scores)
            hits_str = locus.format_hits()
            parent_attrib = (
                f"ID={locus_id};"
                f"Name={locus_id};"
                f"n_hits={locus.n_hits};"
                f"evalues={evalues_str};"
                f"scores={scores_str};"
                f"hits={hits_str};"
                f"size={locus.size};"
            )

            gff3.write(
                "\t".join([
                    locus.target,
                    source,
                    "putative_locus",
                    str(locus.start),
                    str(locus.end),
                    ".",
                    locus.strand,
                    ".",
                    parent_attrib,
                ]) + "\n"
            )

            if span is not None:
                gff3.write(
                    "\t".join([
                        locus.target,
                        source,
                        "putative_locus_span",
                        str(locus.span_start(span)),
                        str(locus.span_end(span)),
                        ".",
                        locus.strand,
                        ".",
                        parent_attrib,
                    ]) + "\n"
                )

            if expansion is not None:
                gff3.write(
                    "\t".join([
                        locus.target,
                        source,
                        "putative_locus_expanded",
                        str(locus.extstart),
                        str(locus.extend),
                        ".",
                        locus.strand,
                        ".",
                        parent_attrib,
                    ]) + "\n"
                )


            for hit_num, hit in enumerate(locus.transcript_ordered_hits, start=1):
                hit_id = f"{locus_id}.hit{hit_num}"
                hit_attrib = (
                    f"ID={hit_id};"
                    f"Parent={locus_id};"
                    f"hmmfrom={hit.hmmfrom};"
                    f"hmmto={hit.hmmto};"
                    f"evalue={hit.evalue:.3g};"
                    f"score={hit.score:.1f}"
                )
                gff3.write(
                    "\t".join([
                        locus.target,
                        source,
                        "putative_exon",
                        str(hit.gstart),
                        str(hit.gend),
                        f"{hit.score:.1f}",
                        locus.strand,
                        ".",
                        hit_attrib,
                    ]) + "\n"
                )


def main():
    """
    Arguments to be passed to the command line.
    """
    p = argparse.ArgumentParser(
        description="Reconstruct putative loci from an ordered nhmmer --tblout output and export to gff3 file."
    )

    p.add_argument(
        "-i", "--input", required=True, help="Ordered nhmmer --tblout file without description column.",
    )

    p.add_argument(
        "-o", "--output", default="putative_loci.gff3", help="gff3 output file for putative loci reconstructed (default: putative_loci.gff3).",
    )

    p.add_argument(
        "-c", "--coord-source", choices=["ali", "env"], default="ali", help="Genomic coordinates to use: alifrom/alito (default: core alignment, narrower) or envfrom/envto (envelope, wider).",
    )

    p.add_argument(
        "-e", "--expansion", type=int, default=None, help="Number of bases to expand around coordinates of putative loci reconstructed.",
    )

    p.add_argument(
        "-s","--span", type=int, default=None, help="Query span in the genome, e.g. expected approximate span of a fully reconstructed loci.",
    )

    p.add_argument(
        "-m","--max-size", type=int, default=None, help="Maximum size (bp) allowed for a reconstructed locus.",
    )

    args = p.parse_args()

    hits = parse_nhmmer_tblout(args.input, coord_source=args.coord_source)
    print(f"nhmmer\n----------\nHits (putative exons) read in {args.input}: {len(hits)}",
        file=sys.stdout)

    merged = merge_overlapping_hits(hits)
    n_merged = len(hits) - len(merged)
    print(
        f"Overlapping hits fused -> {n_merged}"
        f"\nHits remaining: {len(merged)}",
        file=sys.stdout,
    )

    loci = build_loci(merged, max_size=args.max_size)
    print(f"Loci reconstructed: {len(loci)}",
        file=sys.stdout)

    exp_loci = conditional_expansion(loci, expansion=args.expansion)
    print(f"Loci expanded if expansion parameter on.",
        file=sys.stdout)
    
    write_loci_gff3(exp_loci, args.output, span=args.span, expansion=args.expansion)
    print(f"Putative loci written in {args.output}",
        file=sys.stdout)


if __name__ == "__main__":
    main()


