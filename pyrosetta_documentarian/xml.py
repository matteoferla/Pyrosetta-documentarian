import pyrosetta, os
from collections import defaultdict
from typing import List, Optional
from warnings import warn

from .base import BaseDocumentarian

class XMLDocumentarian(BaseDocumentarian):
    test_folder = 'main/tests/integration/tests'
    _xmlfilenames = []
    _xmlprotocols = {}
    _mover_directory = defaultdict(list)

    def get_relevant_scripts(self) -> List[str]:
        """
        Get the paths of the XML test scripts that use the given mover.

        :return:
        """
        mover_name = self.target.get_name()
        if mover_name in self.mover_directory:
            return self.mover_directory[mover_name]
        else:
            return []

    def get_mover_from_script(self, xmlfilename: str, target: Optional=None):
        if target is None:
            target = self.target
        protocol = self.load_xmlfilename(xmlfilename)
        movers = self.get_movers_in_protocol(protocol)
        for mover in movers:
            if mover.get_name() == target.get_name():
                return mover
        else:
            raise ValueError(f'Could not find mover {target.get_name()} in {xmlfilename}')

    @property  # cached
    def mover_directory(self):
        """
        This is a dictionary associates a mover name with a list of files that use it.
        This can be used with xmlprotocols associates a file with a protocol obj
        and then self.get_movers_in_protocol(proto).
        The latter can be compared with the attribute comparison methods.
        """
        if self._mover_directory:
            return self._mover_directory
        for filename, proto in self.xmlprotocols.items():
            for mover in self.get_movers_in_protocol(proto):
                self._mover_directory[mover.get_name()].append(filename)
        return self._mover_directory

    @property  # cached
    def xmlfilenames(self):
        if self._xmlfilenames != []:
            return self._xmlfilenames

        # walk
        def get_files(path, files: list):
            if os.path.isdir(path):
                for file in os.listdir(path):
                    get_files(os.path.join(path, file), files)
            elif '.xml' in path:
                files.append(path)
            else:
                pass
            return files

        self._xmlfilenames = get_files(os.path.join(self.rosetta_folder, self.test_folder), [])
        return self._xmlfilenames

    @property  # cached
    def xmlprotocols(self):
        if self._xmlprotocols != {}:
            return self._xmlprotocols
        for xmlfilename in self.xmlfilenames:
            try:
                self._xmlprotocols[xmlfilename] = self.load_xmlfilename(xmlfilename)
            except Exception as error:
                warn(f'{xmlfilename} could not be read due to {error.__class__.__name__}: {error}')
        return self._xmlprotocols

    def load_xmlfilename(self, xmlfilename: str):
        pose = pyrosetta.Pose()
        xml_obj = pyrosetta.rosetta.protocols.rosetta_scripts.RosettaScriptsParser()
        protocol = xml_obj.generate_mover_and_apply_to_pose(pose, xmlfilename)
        return protocol

    def get_movers_in_protocol(self, protocol: pyrosetta.rosetta.protocols.rosetta_scripts.ParsedProtocol):
        # I could not figure out how to find out how many movers are in a protocol.
        # Namely how many Add tags in PROTOCOLS in the XML
        i = 1
        movers = []
        while True:
            try:
                movers.append(protocol.get_mover(i))
            except RuntimeError:
                break
        # movers may have movers themselves. All the movers in MOVERS
        for mover in movers:
            if hasattr(mover, 'mover'):
                movers.append(mover.mover())
        return mover

