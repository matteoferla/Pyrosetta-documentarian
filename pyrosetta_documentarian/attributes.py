import pyrosetta
import pandas as pd
from typing import Tuple, List, Dict, Set, Any, Optional, Sequence

from .base import BaseDocumentarian

class AttributeDocumentarian(BaseDocumentarian):
    """
    Analyses a Pyrosetta object and determines what is different from default.
    For example. Give a working XML script:

    >>> xml_obj = pyrosetta.rosetta.protocols.rosetta_scripts.RosettaScriptsParser()
    >>> protocol = xml_obj.generate_mover_and_apply_to_pose(pose, 'script.xml')
    >>> protocol.apply(pose)

    One can reverse engineer it, thusly:

    >>> pm = protocol.get_mover(1)
    >>> print(pm.mover_name()) # mover called in script!
    >>> AttributeDocumentarian(pm).compare(evo) # -> pd.DataFrame

    ---------------------------
    Attributes:

    * target: instance
    * target_cls: class
    * base: The tuple of classes inherited (``__mro__``)
    * uninherited: The set of attributes that are absent in the parent class
    * citation: string of citation

    Methods:

    * describe(): describe attributes
    * test(): calls the methods
    * compare(): compares the results of a ``test()`` to that of a blank instance


    """

    @property
    def uninherited(self) -> Set[str]:
        """
        The set of attributes that are absent in the parent class.
        Has no idea if other were overwritten though!

        :rtype: Set[str]
        """
        if len(self.base) > 1:
            return set(dir(self.base[0])) - set(dir(self.base[1]))

    def describe(self, iterable: Optional[Sequence[str]] = None) -> None:
        """
        Describe attributes by calling help.
        If ``iterable`` is provided, it will print only those.
        """
        if iterable is None:
            iterable = dir(self.target)
        for methodname in iterable:
            print(f'## {methodname}')
            method = getattr(self.target, methodname)
            help(method)

    def test(self,
             iterable: Optional[Sequence[str]] = None,
             silent: bool = True) -> Dict[str, Any]:
        """
        Calls without arguments the methods.
        If ``iterable`` is provided, it will call only those.
        Returns a dictionary of the results.
        """
        if iterable is None:
            iterable = dir(self.target)
        results = {}
        for methodname in iterable:
            method = getattr(self.target, methodname)
            try:
                result = method()
                results[methodname] = result
                if silent is False:
                    print(f'Calling worked for {methodname}: {result}')
            except TypeError as error:
                results[methodname] = 'N/A'
                if silent is False:
                    print(f'Calling failed for {methodname}: {result}')
        return results

    def test_uninherited(self, silent: bool = True) -> dict:
        """
        Calls without arguments the methods that where not inherited.
        """
        return self.test(self.uninherited, silent)

    def compare(self, reference: Optional[pyrosetta.rosetta.protocols.moves.Mover] = None) -> pd.DataFrame:
        """
        Tests the methods (see ``test()`` and compares them to a generic instance
        or to ``reference`` if provided.
        """
        c = self.test()
        if reference is None:
            reference = self.target_cls()
        refexplorer = self.__class__(reference)
        r = refexplorer.test()
        return self._make_table(c, r)

    def compare_uninherited(self, reference: Optional[pyrosetta.rosetta.protocols.moves.Mover] = None) -> pd.DataFrame:
        """
        Tests the uninherited methods (see ``test()`` and compares them to a generic instance
        or to ``reference`` if provided.
        """
        c = self.test_uninherited()
        if reference is None:
            reference = self.target_cls()
        refexplorer = self.__class__(reference)
        r = refexplorer.test_uninherited()
        return self._make_table(c, r)

    def _make_table(self, case: Dict[str, Any], ref: Dict[str, Any]) -> pd.DataFrame:
        assert case, f'make_table cannot make a table without data (case={case}, ref={ref})'
        proto = [{'attribute': k,
                  'target': ref[k],
                  'reference':  case[k],
                  'different': str(ref[k]) == str(case[k])} for k in case.keys()]
        comparison = pd.DataFrame(proto)
        return comparison.set_index(['attribute'])
