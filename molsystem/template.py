# -*- coding: utf-8 -*-
import logging

from .align import AlignMixin
from .cif import CIFMixin
from .cms_schema import CMSSchemaMixin
from .inchi import InChIMixin
from .molfile import MolFileMixin
from .openbabel import OpenBabelMixin
from .pubchem import PubChemMixin
from .rdkit_ import RDKitMixin
from .pdb import PDBMixin
from .qcschema import QCSchemaMixin
from .smiles import SMILESMixin
from .topology import TopologyMixin

logger = logging.getLogger(__name__)


class _Template(
    PDBMixin,
    MolFileMixin,
    CIFMixin,
    CMSSchemaMixin,
    InChIMixin,
    SMILESMixin,
    TopologyMixin,
    OpenBabelMixin,
    PubChemMixin,
    RDKitMixin,
    QCSchemaMixin,
    AlignMixin,
    object,
):
    """A class providing the API for templates.

    There are two types of templates:

    simple
        Which have a category and name, which together are unique.
    full
        Which have a category and name, which together are unique,
        and also a reference to a configuration that
        describes the atoms, bonds, and other properties of the
        template.

    Parameters
    ----------
    system_db : SystemDB
        The system database that this template belongs to.
    tid : int
        The id of the template in the template table.
    logger : logging.Logger = logger
        A logger to use in place of the one from this module.
    """

    def __init__(self, system_db, tid, logger=logger):
        self._system_db = system_db
        self._id = tid
        self._logger = logger

        # Cache for performance:
        self._is_full = None
        self._configuration = None
        self._configuration_id = None

    @property
    def category(self):
        """The category of this template."""
        sql = "SELECT category FROM template WHERE id = ?"
        self.cursor.execute(sql, (self.id,))
        return self.cursor.fetchone()[0]

    @property
    def configuration(self):
        """The configuration for a full template.

        Raises
        ------
        TypeError
            If not a full template, so there is no template configuration.

        Returns
        -------
        _Configuration
            The Configuration object.
        """
        if not self.is_full:
            raise TypeError("Not a full template")

        if self._configuration is None:
            cid = self.configuration_id
            if cid is None:
                return None
            self._configuration = self.system_db.get_configuration(cid)
        return self._configuration

    @property
    def configuration_id(self):
        """The configuration id for a full template.

        Raises
        ------
        TypeError
            If not a full template.

        Returns
        -------
        int
            The id of the configuration in the configuration table.
        """
        if not self.is_full:
            raise TypeError("Not a full template")

        if self._configuration_id is None:
            sql = "SELECT configuration FROM template WHERE id = ?"
            self.cursor.execute(sql, (self.id,))
            self._configuration_id = self.cursor.fetchone()[0]
        return self._configuration_id

    @property
    def atoms(self):
        """The atoms for this template.

        Raises
        ------
        TypeError
            If not a full template.

        Returns
        -------
        _Atoms
            The atoms for the template configuration.
        """

        return self.configuration.atoms

    @property
    def bonds(self):
        """The bonds for this template.

        Raises
        ------
        TypeError
            If not a full template.

        Returns
        -------
        _Bonds
            The bonds for the template configuration.
        """

        return self.configuration.bonds

    @property
    def cell(self):
        """The cell for this template.

        Raises
        ------
        TypeError
            If not a full template or is not periodic.

        Returns
        -------
        _Cell
            The cell for the template configuration.
        """

        return self.configuration.cell

    @property
    def coordinate_system(self):
        """The coordinate system for this template.

        Raises
        ------
        TypeError
            If not a full template.

        Returns
        -------
        str
            The coordinate system for the template configuration.
        """

        return self.configuration.coordinate_system

    @property
    def cursor(self):
        """The a cursor for the database."""
        return self.system_db.cursor

    @property
    def db(self):
        """The database connection."""
        return self.system_db.db

    @property
    def density(self):
        """The density for this template.

        Raises
        ------
        TypeError
            If not a full template or is not periodic.

        Returns
        -------
        float
            The density for the template configuration.
        """

        return self.configuration.density

    @property
    def formula(self):
        """The chemical formula for this template.

        Raises
        ------
        TypeError
            If not a full template.

        Returns
        -------
        tuple(str, str, int)
            The chemical formula, empirical formula and Z.
        """

        return self.configuration.formula

    @property
    def is_full(self):
        """Whether this has a template configuration, i.e. is a full template."""
        if self._is_full is None:
            if self._configuration_id is None:
                sql = "SELECT configuration FROM template WHERE id = ?"
                self.cursor.execute(sql, (self.id,))
                cid = self.cursor.fetchone()[0]
                self._is_full = cid is not None
            else:
                self._is_full = True
        return self._is_full

    @property
    def id(self):
        """The id for this template in the template table."""
        return self._id

    @property
    def mass(self):
        """The atomic mass of this template.

        Raises
        ------
        TypeError
            If not a full template.

        Returns
        -------
        float
            The atomic mass of the template configuration.
        """

        return self.configuration.mass

    @property
    def name(self):
        """The name of this template."""
        sql = "SELECT name FROM template WHERE id = ?"
        self.cursor.execute(sql, (self.id,))
        return self.cursor.fetchone()[0]

    @property
    def n_atoms(self):
        """The number of atoms for this template.

        Raises
        ------
        TypeError
            If not a full template.

        Returns
        -------
        int
            The number of atoms for the template configuration.
        """

        return self.configuration.n_atoms

    @property
    def n_bonds(self):
        """The number of bonds for this template.

        Raises
        ------
        TypeError
            If not a full template.

        Returns
        -------
        int
            The number of bonds for the template configuration.
        """

        return self.configuration.n_bonds

    @property
    def periodicity(self):
        """The periodicity of this template.

        Raises
        ------
        TypeError
            If not a full template.

        Returns
        -------
        int
            The periodicity of the template configuration.
        """

        return self.configuration.periodicity

    @property
    def symmetry(self):
        """The symmetry of this template.

        Raises
        ------
        TypeError
            If not a full template.

        Returns
        -------
        _Symmetry
            The symmetry of the template configuration.
        """

        return self.configuration.symmetry

    @property
    def system_db(self):
        """The system_db that we belong to."""
        return self._system_db

    @property
    def volume(self):
        """The volume of this template.

        Raises
        ------
        TypeError
            If not a full template or not periodic.

        Returns
        -------
        float
            The volume of the template configuration.
        """

        return self.configuration.volume
