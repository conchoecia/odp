"""Tests for source/odp_color_manager.py.

Focus on the pure-function color utilities — hex↔rgb conversion, hex
validation, contrasting-color picker, brightness / saturation math,
random / palette color generators. These run without any I/O, so they
exercise the helpers used everywhere downstream (synteny plots, ribbon
plots, ALG colorings).
"""
from __future__ import annotations

import pytest


@pytest.fixture(scope="module")
def cm(source_dir):
    """Import source/odp_color_manager.py via the source_dir fixture."""
    import odp_color_manager
    return odp_color_manager


# ---------------------------------------------------------------------------
# h2r — hex string → normalized RGB triple (each in [0,1])
# ---------------------------------------------------------------------------


def test_h2r_black(cm):
    r, g, b = cm.h2r("#000000")
    assert (r, g, b) == (0.0, 0.0, 0.0)


def test_h2r_white(cm):
    r, g, b = cm.h2r("#FFFFFF")
    assert r == pytest.approx(1.0)
    assert g == pytest.approx(1.0)
    assert b == pytest.approx(1.0)


def test_h2r_red(cm):
    r, g, b = cm.h2r("#FF0000")
    assert r == pytest.approx(1.0)
    assert g == 0.0
    assert b == 0.0


def test_h2r_handles_lowercase(cm):
    assert cm.h2r("#ff0000") == cm.h2r("#FF0000")


def test_h2r_handles_no_hash_prefix(cm):
    assert cm.h2r("ff0000") == cm.h2r("#FF0000")


def test_h2r_known_value(cm):
    # 0x33aa66 -> (51, 170, 102) -> /255 each
    r, g, b = cm.h2r("#33aa66")
    assert r == pytest.approx(51 / 255)
    assert g == pytest.approx(170 / 255)
    assert b == pytest.approx(102 / 255)


# ---------------------------------------------------------------------------
# inverse_color — given a background hex, returns "#000000" or "#ffffff"
# ---------------------------------------------------------------------------


def test_inverse_color_black_returns_black(cm):
    """Black-on-black is the documented short-circuit; lets callers know
    no text was rendered."""
    assert cm.inverse_color("#000000") == "#000000"


def test_inverse_color_white_returns_black(cm):
    """Bright background → dark text."""
    assert cm.inverse_color("#FFFFFF") == "#000000"


def test_inverse_color_dark_returns_white(cm):
    """Dark background → light text."""
    assert cm.inverse_color("#222222") == "#ffffff"


def test_inverse_color_red_returns_white(cm):
    # red has weighted brightness 255*0.299 ≈ 76 → < 186 → white
    assert cm.inverse_color("#FF0000") == "#ffffff"


def test_inverse_color_yellow_returns_black(cm):
    # yellow has weighted brightness ≈ 226 → > 186 → black
    assert cm.inverse_color("#FFFF00") == "#000000"


# ---------------------------------------------------------------------------
# is_valid_hex_code — strict 7-char `#XXXXXX` form
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("good", [
    "#000000", "#FFFFFF", "#abcdef", "#012345", "#ABCDEF",
])
def test_is_valid_hex_code_accepts_well_formed(cm, good):
    assert cm.is_valid_hex_code(good) is True


@pytest.mark.parametrize("bad", [
    "",            # empty
    "#",           # just hash
    "000000",      # missing hash
    "#00000",      # 6 chars total — too short
    "#0000000",    # 8 chars total — too long
    "#GGGGGG",     # non-hex chars
    "#00 000",     # whitespace inside
    "  #000000",   # leading whitespace
    "#000000 ",    # trailing whitespace
    None,          # not a string at all
])
def test_is_valid_hex_code_rejects_malformed(cm, bad):
    if bad is None:
        # Function uses len(); None has no len → TypeError. Accept that
        # too as "rejection".
        with pytest.raises(TypeError):
            cm.is_valid_hex_code(bad)
    else:
        assert cm.is_valid_hex_code(bad) is False


# ---------------------------------------------------------------------------
# calculate_brightness — int color → weighted brightness
# ---------------------------------------------------------------------------


def test_calculate_brightness_black_is_zero(cm):
    assert cm.calculate_brightness(0x000000) == 0


def test_calculate_brightness_white_at_max(cm):
    # 255 * (0.299 + 0.587 + 0.114) = 255 → weighted sum / 1000 = ...
    # Actually formula: (r*299 + g*587 + b*114) / 1000
    expected = (255 * 299 + 255 * 587 + 255 * 114) / 1000
    assert cm.calculate_brightness(0xFFFFFF) == pytest.approx(expected)


def test_calculate_brightness_red(cm):
    expected = (255 * 299) / 1000
    assert cm.calculate_brightness(0xFF0000) == pytest.approx(expected)


def test_calculate_brightness_returns_positive(cm):
    """Any non-zero color should give a positive brightness."""
    assert cm.calculate_brightness(0x010101) > 0


def test_calculate_brightness_is_monotone_in_each_channel(cm):
    """Raising any single channel must not decrease brightness."""
    base = cm.calculate_brightness(0x808080)
    redder = cm.calculate_brightness(0xFF8080)
    greener = cm.calculate_brightness(0x80FF80)
    bluer = cm.calculate_brightness(0x8080FF)
    assert redder >= base
    assert greener >= base
    assert bluer >= base


# ---------------------------------------------------------------------------
# calculate_saturation — int color → [0, 1]
# ---------------------------------------------------------------------------


def test_calculate_saturation_grayscale_is_zero(cm):
    """Gray colors have zero saturation (max == min)."""
    assert cm.calculate_saturation(0x000000) == 0
    assert cm.calculate_saturation(0x808080) == 0
    assert cm.calculate_saturation(0xFFFFFF) == 0


def test_calculate_saturation_pure_red_is_one(cm):
    """Pure red: max=255, min=0 → saturation = 1."""
    assert cm.calculate_saturation(0xFF0000) == pytest.approx(1.0)


def test_calculate_saturation_in_unit_range(cm):
    for c in [0x123456, 0xABCDEF, 0x010101, 0xFEDCBA]:
        s = cm.calculate_saturation(c)
        assert 0.0 <= s <= 1.0


def test_calculate_saturation_black_handles_zero_max(cm):
    """The implementation special-cases max==0 to avoid /0."""
    assert cm.calculate_saturation(0x000000) == 0


# ---------------------------------------------------------------------------
# return_random_color — hex string with leading "#"
# ---------------------------------------------------------------------------


def test_return_random_color_format(cm):
    import re
    for _ in range(20):
        c = cm.return_random_color()
        assert isinstance(c, str)
        assert re.match(r"^#[0-9a-f]{6}$", c), f"bad format: {c!r}"


# ---------------------------------------------------------------------------
# generate_random_color — bounded brightness + saturation
# ---------------------------------------------------------------------------


def test_generate_random_color_is_within_bounds(cm):
    """Sample 30 colors and verify each falls inside the brightness and
    saturation windows the function advertises."""
    for _ in range(30):
        hex_c = cm.generate_random_color()
        assert hex_c.startswith("#") and len(hex_c) == 7
        n = int(hex_c.lstrip("#"), 16)
        b = cm.calculate_brightness(n)
        s = cm.calculate_saturation(n)
        assert 127 < b < 300
        assert 0.3 < s < 0.8


# ---------------------------------------------------------------------------
# yield_color — deterministic palette generator
# ---------------------------------------------------------------------------


def test_yield_color_mypals25_cycles_through_all_25(cm):
    gen = cm.yield_color("mypals25")
    palette = [next(gen) for _ in range(25)]
    assert len(set(palette)) == 25  # all distinct on first cycle
    # 26th value loops back to first
    assert next(gen) == palette[0]


def test_yield_color_tableau20_default_length(cm):
    gen = cm.yield_color("tableau20")
    palette = [next(gen) for _ in range(20)]
    assert len(set(palette)) == 20
    assert next(gen) == palette[0]


def test_yield_color_unknown_palette_raises(cm):
    with pytest.raises(KeyError):
        # KeyError comes from the dict lookup before the generator yields.
        next(cm.yield_color("does-not-exist"))


def test_yield_color_each_entry_is_valid_hex(cm):
    gen = cm.yield_color("mypals25")
    import re
    for _ in range(25):
        c = next(gen)
        assert re.match(r"^#[0-9A-Fa-f]{6}$", c), f"bad palette entry: {c!r}"
