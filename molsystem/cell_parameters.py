# -*- coding: utf-8 -*-

import logging

from molsystem.cell import Cell
from molsystem.table import _Table as Table

logger = logging.getLogger(__name__)

labels = ('a', 'b', 'c', 'alpha', 'beta', 'gamma')


class _CellParameters(Table):
    """The representation of the periodic cell
    """

    def __init__(self, system, table='cell'):

        super().__init__(system, table)

        self._configuration_table = self._system['configuration']

    def cell(self, configuration=None):
        """Return the cell parameters for the configuration.

        Parameters
        ----------
        configuration : int = None
            The configuration of interest. Defaults to the current
            configuration.

        Returns
        -------
        cell : Cell
            The cell parameters as a Cell object.
        """
        cell_id, configuration = self.cell_id(configuration)
        if cell_id is None:
            return None
        else:
            self.cursor.execute(
                f"SELECT a, b, c, alpha, beta, gamma FROM {self.table}"
                "  WHERE id = ?", (cell_id,)
            )
            return Cell(*self.cursor.fetchone())

    def cell_id(self, configuration=None):
        """Return the id of the cell parameters for the configuration.

        Parameters
        ----------
        configuration : int = None
            The configuration of interest. Defaults to the current
            configuration.

        Returns
        -------
        cell_id : int
            The id of the cell parameters.
        """
        if configuration is None:
            configuration = self.system.current_configuration
        self.cursor.execute(
            "SELECT cell FROM configuration WHERE id = ?", (configuration,)
        )
        cell_id = self.cursor.fetchone()
        if cell_id is None:
            return None, configuration
        else:
            return cell_id[0], configuration

    def set_cell(self, *args, configuration=None):
        """Set the cell parameters for the configuration.

        Parameters
        ----------
        args : Cell, iterable or 6 floats
            The cell parameters as a Cell object, an iterable of length 6
            (list, tuple, ...), or 6 floats

        configuration : int = None
            The configuration of interest. Defaults to the current
            configuration.

        Returns
        -------
        None
        """
        if len(args) == 1:
            if isinstance(args[0], Cell):
                parameters = args[0].parameters
            elif len(args[0]) == 6:
                parameters = [*args[0]]
            else:
                raise ValueError(
                    'cell must be a 6-vector or six separate values'
                )
        elif len(args) == 6:
            parameters = [*args]
        else:
            raise ValueError('cell must be a 6-vector or six separate values')

        cell_id, configuration = self.cell_id(configuration)
        if cell_id is None:
            a, b, c, alpha, beta, gamma = parameters
            cell_id = self.append(
                a=a, b=b, c=c, alpha=alpha, beta=beta, gamma=gamma
            )[0]
            self.cursor.execute(
                "UPDATE configuration SET cell = ? WHERE id = ?",
                (cell_id, configuration)
            )
        else:
            parameters.append(cell_id)
            self.cursor.execute(
                "UPDATE cell SET a=?, b=?, c=?, alpha=?, beta=?, gamma=?"
                " WHERE id = ?", parameters
            )
