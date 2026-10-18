"""Tabulated data about the elements."""

# The element data, with the atomic weights, are seamm_util's, which takes them from
# mendeleev; this module keeps the lookups that molsystem and the steps use.
from seamm_util.elemental_data import element_data as data

symbol_to_atno = {}
atno_to_symbol = {}
symbol_to_mass = {}
atno_to_mass = {}

for symbol, value in data.items():
    symbol_to_atno[symbol] = value["atomic number"]
    atno_to_symbol[value["atomic number"]] = symbol
    symbol_to_mass[symbol] = value["atomic weight"]
    atno_to_mass[value["atomic number"]] = value["atomic weight"]

del symbol
del value


def to_atnos(symbols):
    """Convert element symbols to atomic numbers.

    Parameters
    ----------
    symbols : [str]
        The atomic symbols

    Returns
    -------
    atnos : [int]
        The corresponding atomic numbers (1..118)
    """
    return [symbol_to_atno[x] for x in symbols]


def to_symbols(atnos):
    """Convert atomic numbers to element symbols.

    Parameters
    ----------
    atnos : [int]
        The atomic numbers (1..118)

    Returns
    -------
    symbols : [str]
        The corresponding atomic symbols
    """
    return [atno_to_symbol[x] for x in atnos]


def masses(atno_or_symbols):
    """Get the atomic mass given atomic symbols or numbers.

    Parameters
    ----------
    atno_or_symbols : [int or str]
        The atomic numbers or symbols

    Returns
    -------
    masses : [float]
        The atomic masses
    """
    result = []
    for value in atno_or_symbols:
        if isinstance(value, int):
            result.append(atno_to_mass[value])
        else:
            result.append(symbol_to_mass[value])

    return result
