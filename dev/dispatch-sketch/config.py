"""
Shared coverage configuration.

The whole point of the dispatcher is that BOTH backends (the hand-rolled pysam
calculator and mosdepth) consume this one object. Anything backend-specific
lives behind the backend boundary, not here. Where a knob can't be honored by a
backend, that backend should raise -- silent divergence between backends is the
worst outcome, since the main reason to have both is cross-validation.

Flag masks mirror samtools/mosdepth integer semantics. Note the default-mask
mismatch between tools, spelled out below: we pin an explicit shared default so
the two backends actually agree out of the box.
"""

"""
mosdepth_config:
  include_contigs:
    - "chr*"
  exclude_contigs:
    - "*_alt"
    - "*_decoy"
    - "*_random"
    - "chrUn*"
    - "HLA*"
    - "chrM"
    - "chrEBV"
"""

from dataclasses import dataclass
from pathlib import Path

# SAM flag name -> bit (kept here so config/CLI can speak names, internals use ints).
SAM_FLAGS = {
    "paired":        0x1,
    "proper_pair":   0x2,
    "unmapped":      0x4,
    "mate_unmapped": 0x8,
    "reverse":       0x10,
    "mate_reverse":  0x20,
    "read1":         0x40,
    "read2":         0x80,
    "secondary":     0x100,
    "qcfail":        0x200,
    "duplicate":     0x400,
    "supplementary": 0x800,
}

# Explicit shared default = unmapped | secondary | qcfail | duplicate | supplementary.
# This is deliberately the UNION of the two tools' native defaults so neither
# backend is silently doing something different:
#   - read_filter.py default   = unmapped|secondary|duplicate|supplementary  (3332)
#   - mosdepth --flag default   = unmapped|secondary|qcfail|duplicate        (1796)
# Pinning it here (3844) removes that footgun; override per-run if you want to
# reproduce a specific tool's out-of-the-box number.
DEFAULT_EXCLUDE = (
    SAM_FLAGS["unmapped"]
    | SAM_FLAGS["secondary"]
    | SAM_FLAGS["qcfail"]
    | SAM_FLAGS["duplicate"]
    | SAM_FLAGS["supplementary"]
)


def to_mask(flags) -> int:
    """Normalize a flag spec (int mask | iterable of names/ints | None) to an int."""
    if flags is None:
        return 0
    if isinstance(flags, int):
        return flags
    mask = 0
    for f in flags:
        mask |= f if isinstance(f, int) else SAM_FLAGS[f]
    return mask


@dataclass
class CoverageConfig:
    # --- inputs (shared) ---
    cram: Path
    reference: Path
    index: Path | None = None            # defaults to <cram>.crai if None

    # --- backend selection ---
    backend: str = "native"              # "native" | "mosdepth"

    # Optional explicit contig subset (None = all). native loops only these;
    # mosdepth would need a --chrom (single) or a BED to match.
    chroms: tuple[str, ...] | None = None

    # Glob include/exclude over reference names, applied when `chroms` is None.
    # Keeps decoy/unplaced artifact contigs out of the headline report.
    include_contigs: tuple[str, ...] = ("chr*",)
    exclude_contigs: tuple[str, ...] = ("*_alt", "*_random", "chrUn*", "*_decoy", "HLA*", "chrEBV")

    # --- read filtering (shared knobs) ---
    min_mapping_quality: int = 0         # native -Q  /  mosdepth --mapq
    include_flags: int = 0               # native -f  /  mosdepth --include-flag
    exclude_flags: int = DEFAULT_EXCLUDE  # native -F  /  mosdepth --flag
    exclude_all_flags: int = 0           # native -G  /  mosdepth: NO EQUIVALENT (backend raises)

    # --- outputs / execution ---
    per_base: bool = False               # keep full per-base depth vector (native) /
                                         #   drop mosdepth --no-per-base when True
    threads: int = 1                     # mosdepth --threads; native is single-threaded

    def __post_init__(self):
        self.cram = Path(self.cram)
        self.reference = Path(self.reference)
        self.index = Path(self.index) if self.index else self.cram.with_suffix(self.cram.suffix + ".crai")
        self.include_flags = to_mask(self.include_flags)
        self.exclude_flags = to_mask(self.exclude_flags)
        self.exclude_all_flags = to_mask(self.exclude_all_flags)

    def select_contigs(self, references) -> list[str]:
        """Resolve which contigs to process: explicit `chroms` wins, else apply
        include-then-exclude globs over the CRAM's reference names."""
        import fnmatch
        if self.chroms is not None:
            return list(self.chroms)
        kept = [r for r in references
                if any(fnmatch.fnmatch(r, p) for p in self.include_contigs)]
        return [r for r in kept
                if not any(fnmatch.fnmatch(r, p) for p in self.exclude_contigs)]
