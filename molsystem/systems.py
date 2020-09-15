# -*- coding: utf-8 -*-

"""A dictionary-like object for holding systems
"""

import collections.abc
import logging
from pathlib import Path
import shutil
import sqlite3
import tempfile

import pathvalidate

from molsystem.system import _System

logger = logging.getLogger(__name__)


class Systems(collections.abc.MutableMapping):

    def __init__(self):
        """Initialize the systems with no systems."""

        self._systems = {}

    def __getitem__(self, key):
        """Allow [] access to the dictionary of systems"""
        return self._systems[key]['system']

    def __setitem__(self, key, value):
        """Allow x[key] access to the data"""
        raise NotImplementedError(f"Table '{key}' cannot be created yet")

    def __delitem__(self, key):
        """Allow deletion of keys"""
        if key in self:
            data = self._systems[key]
            del data['system']
            if data['temporary']:
                if 'tempdir' in data:
                    shutil.rmtree(data['tempdir'])
                else:
                    data['path'].unlink()
            del self._systems[key]
        else:
            raise KeyError(
                f"Trying to delete system '{key}', which does not exist."
            )

    def __iter__(self):
        """Allow iteration over the object"""
        return iter(self._systems)

    def __len__(self):
        """The len() command"""
        return len(self._systems)

    def __repr__(self):
        """The string representation of this object"""
        return repr(self._systems)

    def __str__(self):
        """The pretty string representation of this object"""
        return str(self._systems)

    def __contains__(self, key):
        """Return a boolean indicating if a system exists."""
        return key in self._systems

    def __eq__(self, other):
        """Return a boolean if this object is equal to another"""
        raise NotImplementedError()

    def list(self):
        """Return a list of the systems."""
        return [*self._systems]

    def create_system(self, name, filename=None, temporary=False, force=False):
        """Create a system with a given name, and optionally a filename."""
        if name in self:
            raise KeyError(f"System '{name}' already exists.")

        self._systems[name] = data = {}
        data['temporary'] = temporary

        if temporary:
            tmp = pathvalidate.sanitize_filename(name + '.db', platform='auto')
            data['tempdir'] = Path(tempfile.mkdtemp())
            path = data['tempdir'] / tmp
        elif filename is None:
            tmp = pathvalidate.sanitize_filename(name + '.db', platform='auto')
            path = Path(tmp)
        else:
            tmp = pathvalidate.sanitize_filename(filename, platform='auto')
            path = Path(tmp)

        path = path.expanduser().resolve()
        if path.exists():
            if force:
                path.unlink()
            else:
                raise RuntimeError(f"File '{path}' exists!")

        filename = str(path)
        system = _System(self, nickname=name, filename=filename)

        data['system'] = system
        data['path'] = path

        return system

    def copy_system(
        self, other, name=None, filename=None, temporary=False, force=False
    ):
        """Create a copy of a system, optionally with a given name and
        filename."""

        if name is None:
            tmp_name = other.name
            name = tmp_name + '_copy'
            i = 1
            while name in self:
                i += 1
                name = f'{tmp_name}_copy_{i}'
        elif name in self:
            raise KeyError(f"System '{name}' already exists.")

        self._systems[name] = data = {}
        data['temporary'] = temporary

        if temporary:
            tmp = pathvalidate.sanitize_filename(name + '.db', platform='auto')
            data['tempdir'] = Path(tempfile.mkdtemp())
            path = data['tempdir'] / tmp
        elif filename is None:
            tmp = pathvalidate.sanitize_filename(name + '.db', platform='auto')
            path = Path(tmp)

        # Get simple absolute path
        path = path.expanduser().resolve()

        if path.exists():
            if force:
                path.unlink()
            else:
                raise RuntimeError(f"File '{path}' exists!")

        filename = str(path)

        # Copy the database over
        db = sqlite3.connect(filename)
        other.db.commit()
        other.db.backup(db)
        db.close()

        # and open it
        system = _System(self, nickname=name, filename=filename)

        data['system'] = system
        data['path'] = path

        return system

    def open_system(self, filename, name=None, temporary=False):
        """Open an existing system with a given name."""
        path = Path(filename)
        if name is None:
            name = path.with_suffix('').name

        if name in self:
            raise RuntimeError(f"System '{name}' already exists!")

        path = path.expanduser().resolve()
        if not path.exists():
            raise RuntimeError(f"File '{path}' does not exist!")

        self._systems[name] = data = {}
        data['temporary'] = temporary

        filename = str(path)

        # and open it
        system = _System(self, nickname=name, filename=filename)

        data['system'] = system
        data['path'] = path

        return system

    def overwrite(self, system, other):
        """Overwrite a system with the contents of another."""

        system.db.commit()
        system.cursor.close()
        system.db.close()
        path = self._systems[system.nickname]['path']
        path.unlink()
        system._db = sqlite3.connect(path)
        system._db.row_factory = sqlite3.Row
        system._db.execute('PRAGMA foreign_keys = ON')
        system._cursor = system._db.cursor()
        other.db.commit()
        other.db.backup(system._db)
