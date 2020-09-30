# -*- coding: utf-8 -*-

"""Functions for handling SMILES"""

import logging

from openbabel import openbabel

logger = logging.getLogger(__name__)


class SMILESMixin:
    """A mixin for handling SMILES."""

    def to_smiles(
        self,
        canonical=False,
        configuration=None,
        name=False,
        hydrogens=False
    ):
        """Create the SMILES string from the system.

        Parameters
        ----------
        canonical : bool = False
            Whether to create canonical SMILES
        configuration : int = None
            The configuration to use, defaults to the current configuration.
        name : bool = False
            Whether to return the name of the system
        hydrogens : bool = False
            Whether to keep H's in the SMILES string.

        Returns
        -------
        str
            The SMILES string, or (SMILES, name) if the rname is requested
        """
        logger.info('to_smiles')

        molfile = self.to_molfile_text()
        logger.info('...molfile = ')
        logger.info(molfile)
        logger.info('...end molfile')

        obConversion = openbabel.OBConversion()
        if canonical:
            obConversion.SetInAndOutFormats("mdl", "can")
        else:
            obConversion.SetInAndOutFormats("mdl", "smi")
        mol = openbabel.OBMol()
        obConversion.ReadString(mol, molfile)

        logger.info(f'MolFile has {mol.NumAtoms()} atoms')

        if not name:
            obConversion.AddOption('n')
        if hydrogens:
            obConversion.AddOption('h')
        smiles = obConversion.WriteString(mol)

        logger.info('')
        logger.info('')
        logger.info('')
        logger.info('')
        logger.info(f"smiles = '{smiles}'")

        if name:
            return smiles.strip().split('\t')
        else:
            return smiles.strip()

    def from_smiles(self, smiles, name=None, configuration=None):
        """Create the system from a SMILES string.

        Parameters
        ----------
        smiles : str
            The SMILES string
        name : str = None
            The name of the molecule
        configuration : int = None
            The configuration to use, defaults to the current configuration.

        Returns
        -------
        None
        """

        obConversion = openbabel.OBConversion()
        obConversion.SetInAndOutFormats("smi", "mdl")
        obConversion.AddOption('3')
        mol = openbabel.OBMol()
        obConversion.ReadString(mol, smiles)

        # Add hydrogens
        mol.AddHydrogens()

        # Get coordinates for a 3-D structure
        builder = openbabel.OBBuilder()
        builder.Build(mol)

        molfile = obConversion.WriteString(mol)
        self.from_molfile_text(molfile)
        if name is not None:
            self.name = name
