"""A deterministic PRNG that is bit-identical in Python and JavaScript.

The two frustration engines in this project are cross-validated against each
other. That check is only worth anything if both draw the *same* decoys, so
the generator cannot be `random` on one side and `Math.random` on the other.

mulberry32 is used because it is a handful of 32-bit integer operations, which
both languages can reproduce exactly: Python masks to 32 bits, JavaScript gets
the same from `|0` and `Math.imul`.
"""

from __future__ import annotations

MASK32 = 0xFFFFFFFF


def _imul(a: int, b: int) -> int:
    """32-bit integer multiply, matching JavaScript's Math.imul."""
    return ((a & MASK32) * (b & MASK32)) & MASK32


class Mulberry32:
    """Seeded 32-bit PRNG. `next_float()` returns a double in [0, 1)."""

    __slots__ = ("_a",)

    def __init__(self, seed: int = 0) -> None:
        self._a = seed & MASK32

    def next_float(self) -> float:
        self._a = (self._a + 0x6D2B79F5) & MASK32
        t = self._a
        t = _imul(t ^ (t >> 15), 1 | t)
        t = (t + _imul(t ^ (t >> 7), 61 | t)) & MASK32 ^ t
        return ((t ^ (t >> 14)) & MASK32) / 4294967296.0

    def next_int(self, n: int) -> int:
        """Uniform integer in [0, n)."""
        return int(self.next_float() * n)
