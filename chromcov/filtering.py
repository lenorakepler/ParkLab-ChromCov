"""
Read filtering: SAM flag vocabulary + the per-read `ReadFilter`.

These are SAM *read* flags (unmapped/duplicate/secondary/...) and the samtools
`-Q/-f/-F/-G` filtering they drive. They are kept OUT of the config model (which
only speaks field values) and separate from the coverage-QC *abnormality* flags
in `policy.py` -- the two senses of "flag" live in different modules on purpose.
"""
from __future__ import annotations

from dataclasses import dataclass

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

# Default read exclusions: unmapped | secondary | qcfail | supplementary  (2820).
# Duplicates are deliberately KEPT. In high-depth cancer data, duplicate-marking
# flags reads that share start/end coordinates -- but at high depth independent
# fragments increasingly collide on coordinates by chance, so removing "duplicates"
# can under-count genuinely amplified regions. Override with -F to reproduce a
# tool's native default (e.g. mosdepth also excludes duplicates). The mosdepth
# cross-check add-on defaults to THIS same mask, so the two still agree out of the box.
DEFAULT_EXCLUDE = (
    SAM_FLAGS["unmapped"]
    | SAM_FLAGS["secondary"]
    | SAM_FLAGS["qcfail"]
    | SAM_FLAGS["supplementary"]
)

DEFAULT_INCLUDE_CONTIGS = ("chr*",)
DEFAULT_EXCLUDE_CONTIGS = ("*_alt", "*_random", "chrUn*", "*_decoy", "HLA*", "chrEBV")


def to_mask(flags) -> int:
    """
    Outputs a bitmask: either integers, flags (will be converted), or None.
    Quicker filter while still allowing for easy specification in config.

    Raises ValueError on unknown flag
    """
    if flags is None:
        return 0
    if isinstance(flags, int):
        return flags
    mask = 0
    for f in flags:
        if isinstance(f, int):
            mask |= f
        elif f in SAM_FLAGS:
            mask |= SAM_FLAGS[f]
        else:
            raise ValueError(
                f"unknown SAM flag {f!r}; choose from {sorted(SAM_FLAGS)} or pass an int mask"
            )
    return mask


@dataclass
class ReadFilter:
    """
    Build a read filter based on config specification. Create once, use as filter for each read.
    """
    # -f: require a read to have ALL these flags
    include_flags: int = 0

    # -F: exclude a read if it has ANY of these flags
    exclude_flags: int = DEFAULT_EXCLUDE

    # -G: exclude a read only if it has ALL these flags
    exclude_all_flags: int = 0

    # -Q: minimum mapping quality a read must have
    min_mapping_quality: int = 0

    # Create the mask
    def __post_init__(self):
        self.include_flags     = to_mask(self.include_flags)
        self.exclude_flags     = to_mask(self.exclude_flags)
        self.exclude_all_flags = to_mask(self.exclude_all_flags)

    @classmethod
    def from_config(cls, cfg) -> "ReadFilter":
        """Build the read filter from a Config -- the single place both the mean and
        the --full path derive their filter (was the duplicated `_read_filter`)."""
        return cls(
            include_flags=cfg.include_flags,
            exclude_flags=cfg.exclude_flags,
            exclude_all_flags=cfg.exclude_all_flags,
            min_mapping_quality=cfg.min_mapping_quality,
        )

    def fails(self, read):
        flag = read.flag
        if flag & self.exclude_flags:
            return True
        if (flag & self.include_flags) != self.include_flags:
            return True
        if self.exclude_all_flags and (flag & self.exclude_all_flags) == self.exclude_all_flags:
            return True
        if read.mapping_quality < self.min_mapping_quality:
            return True
        return False
