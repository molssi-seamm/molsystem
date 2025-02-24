# -*- coding: utf-8 -*-

"""A dictionary-like object for holding templates"""

from collections.abc import Sequence
import logging

from .table import _Table
from .template import _Template

logger = logging.getLogger(__name__)


class _Templates(_Table):
    """Templates -- examples -- for subsets.

    A template serves two purposes. First, it provides a named
    instance for labeling similar subsets. The name of a template has
    two parts, a category and a name, e.g. residue/ala,
    molecule/water, or general/group1. This provides a way to group
    and identify subsets.

    Second, the template may optionally point to a configuration of a
    system, in which case the atoms in the template configuration are
    templates for the atoms in subsets created for the template. This
    provides a way to create new molecules or fragments from a
    template, or to transfer information from a template into a
    configuration. For example, this could be used to transfer
    canonical names and atom categories residue templates onto a
    protein molecule.

    This class handles the templates.
    """

    def __init__(self, system_db, logger=logger):
        """Initialize the instance."""
        super().__init__(system_db, table="template", logger=logger)

    # def __repr__(self):
    #     """The string representation of this object"""
    #     raise NotImplementedError()

    # def __str__(self):
    #     """The pretty string representation of this object"""
    #     raise NotImplementedError()

    @property
    def categories(self):
        """Return the current categories of templates."""
        sql = "SELECT DISTINCT category FROM template ORDER BY category ASC"
        return [x[0] for x in self.db.execute(sql)]

    @property
    def n_templates(self):
        """The number of templates."""
        return self.n_rows

    def create(
        self, name, category="general", configuration=None, smiles=True, **kwargs
    ):
        """Create one templates.

        This is an explicit wrapper of the Table method append, specializing
        it for creating one template.

        Parameters
        ----------
        name : str
            The name of the template/subset. The name-category pair must be
            unique.
        category : str = 'general'
            An optional category of template/subset, e.g. 'general', 'residue',
            'molecule'
        configuration : int or _Configuration = None
            An optional template configuration or its id to use.
        smiles : bool = True
            Create and store the canonical SMILES for full templates.
        kwargs : keyword arguments
            Other existing attributes and values.

        Returns
        -------
        _template
            The created template.
        """
        if configuration is None:
            cid = None
        else:
            if isinstance(configuration, int):
                cid = configuration
                configuration = self.system_db.get_configuration(cid)
            else:
                cid = configuration.id
            canonical_smiles = configuration.to_smiles(canonical=True)
            kwargs["canonical_smiles"] = canonical_smiles

        tid = self.append(category=category, name=name, configuration=cid, **kwargs)[0]
        return _Template(self.system_db, tid)

    def create_many(
        self, name, category="general", configuration=None, smiles=True, **kwargs
    ):
        """Create one or more templates.

        This is an explicit wrapper of the Table method append, specializing
        it for templates.

        Parameters
        ----------
        name : str
            The name of the template/subset. The name-category pair must be
            unique.
        category : str = 'general'
            An optional category of template/subset, e.g. 'general', 'residue',
            'molecule'
        configuration : int or _Configuration = None
            An optional template configuration or its id to use.
        smiles : bool = True
            Create and store the canonical SMILES for full templates.
        kwargs : keyword arguments
            Other existing attributes and values.

        Returns
        -------
        [_Template]
            The created templates.
        """
        if configuration is None:
            cids = None
        else:
            cids = []
            if isinstance(configuration, Sequence):
                kwargs["canonical_smiles"] = []
                if isinstance(configuration, int):
                    cid = configuration
                    configuration = self.system_db.get_configuration(cid)
                else:
                    cid = configuration.id
                canonical_smiles = configuration.to_smiles(canonical=True)
                kwargs["canonical_smiles"].append(canonical_smiles)
                cids.append[cid]
            else:
                if isinstance(configuration, int):
                    cids = configuration
                    configuration = self.system_db.get_configuration(cids)
                else:
                    cids = configuration.id
                canonical_smiles = configuration.to_smiles(canonical=True)
                kwargs["canonical_smiles"] = canonical_smiles

        tids = self.append(category=category, name=name, configuration=cids, **kwargs)
        return [_Template(self.system_db, tid) for tid in tids]

    def exists(self, name, category="general"):
        """Return whether a given template exists.

        Parameters
        ----------
        name : str
            The name of the template/subset. The name-category pair must be
            unique.
        category : str = 'general'
            An optional category of template/subset, e.g. 'general', 'residue',
            'molecule'

        Returns
        -------
        bool
            Whether the template exists.
        """
        sql = "SELECT COUNT(*) FROM template WHERE category = ? AND name = ?"
        self.cursor.execute(sql, (category, name))
        return self.cursor.fetchone()[0]

    def get(self, name, category="general"):
        """Return a template.

        Parameters
        ----------
        name : str
            The id or name of the template/subset. If the id is given, the
            category is ignored. Otherwise the template with the given name
            and category is returned.
        category : str = 'general'
            An optional category of template/subset, e.g. 'general', 'residue',
            'molecule'

        Returns
        -------
        _Template
            The template, or None if it does not exist.
        """
        if isinstance(name, str):
            # A name is given
            tid = self.get_id(name, category=category)
            if tid is None:
                return None
        else:
            # The id given
            tid = name
        return _Template(self.system_db, tid)

    def get_id(self, name, category="general"):
        """Return the id of a template.

        Parameters
        ----------
        name : str
            The name of the template/subset. The name-category pair must be
            unique.
        category : str = 'general'
            An optional category of template/subset, e.g. 'general', 'residue',
            'molecule'

        Returns
        -------
        int
            The id of the template, or None if it does not exist.
        """
        sql = "SELECT id FROM template WHERE category = ? AND name = ?"
        self.cursor.execute(sql, (category, name))
        result = self.cursor.fetchone()
        if result is not None:
            result = result[0]

        return result

    def names(self, category):
        """Return the current names of templates of the given category.

        Parameters
        ----------
        category : str
            The category of template/subset, e.g. 'general', 'residue',
            'molecule'

        Returns
        -------
        [str]
            The names of templates of the given category.
        """
        sql = "SELECT name FROM template WHERE category = ? ORDER BY name ASC"
        return [x[0] for x in self.db.execute(sql, (category,))]
