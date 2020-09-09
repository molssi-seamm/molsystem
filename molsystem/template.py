# -*- coding: utf-8 -*-

"""A dictionary-like object for holding template


"""

import logging
from typing import TypeVar

from molsystem.table import _Table as Table

System_tp = TypeVar("System_tp", "System", None)
Templates_tp = TypeVar("Templates_tp", "_Templates", str, None)

logger = logging.getLogger(__name__)


class _Template(Table):
    """The Templates class works with SQLite tables describing templates

    See the main documentation of SEAMM for a detailed description of the
    database scheme underlying the system and hence these templates. The
    following tables handle templates:

    template -- A simple table of all of the templates, giving a name and type.

    templateatom -- An optional list of template atoms associated with a
        template.

    templatecoordinates -- An optional set of coordinates for the template
        atoms. These are in Cartesian coordinates since the templates are
        molecular, not periodic, in nature.

    templatebond -- Bonds linking template atoms i & j if present.

    Templates can be added ('append') or removed ('delete').
    """

    def __init__(self, system: System_tp, tablename: str = 'template') -> None:
        super().__init__(system, tablename)

        self._current_template = 1

    @property
    def current_template(self):
        """The template that is the default for the moment."""
        return self._current_template

    @current_template.setter
    def current_template(self, value):
        self.cursor.execute(
            f'SELECT COUNT(*) FROM {self.table} WHERE id = ?', (value,)
        )
        if self.cursor.fetchone()[0] == 0:
            raise KeyError(f"Template '{value}' does not exist.")
        self._current_template = value

    def set_current_template(self, name, type_='general'):
        """Set the current template given the name and optionally type."""
        self.cursor.execute(
            f'SELECT id FROM {self.table} WHERE "name" = ? AND "type" = ?',
            (name, type_)
        )
        row = self.cursor.fetchone()
        if row is None:
            raise KeyError(f"There is no template '{name}' of type '{type_}'.")
        self._current_template = row[0]
