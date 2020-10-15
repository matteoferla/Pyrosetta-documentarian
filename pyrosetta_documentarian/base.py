import pyrosetta
from typing import Tuple


class BaseDocumentarian:
    rosetta_folder = 'rosetta'
    src_folder = 'main/source/src'

    def __init__(self, target: pyrosetta.rosetta.protocols.moves.Mover):
        """
        :param target: a pyrosetta instance to explore
        :type target: pyrosetta.rosetta.protocols.moves.Mover
        """
        if hasattr(target, '__mro__'):  # it's a class-ish (cf. pybind11)
            self.target = target()
            self.target_cls = target
        else:  # it's an instance-ish
            self.target = target
            self.target_cls = target.__class__

    @property
    def base(self) -> Tuple[object]:
        """
        The tuple of classes inherited (``__mro__``)

        :rtype: Tuple[object]
        """
        return self.target_cls.__mro__

    @property
    def citation(self) -> str:
        ci = self.target.provide_citation_info()
        if len(ci) == 0:
            return ''
        else:
            cc = ci.pop()
            # print(cc.module_name(), cc.module_type())
            buffer = pyrosetta.rosetta.std.stringbuf()
            cc.get_citations_formatted(pyrosetta.rosetta.std.ostream(buffer))
            return buffer.str().strip()