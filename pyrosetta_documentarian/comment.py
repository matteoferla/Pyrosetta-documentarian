import os, re, pyrosetta
from collections import defaultdict
from typing import Optional

from .base import BaseDocumentarian


class CommentDocumentarian(BaseDocumentarian):

    def __init__(self, target: pyrosetta.rosetta.protocols.moves.Mover):
        """
        Reads the C++ header file and does a crude extraction of comments.
        """
        super().__init__(target)
        self._comments = defaultdict(str)
        self._arguments = defaultdict(str)
        self._access = defaultdict(str)
        # temps
        self._previous_comments = ''
        self._ongoing_open_parenthesis = ''
        self._multiline = False
        self._current_class = None
        self._current_access = 'private'  # private | public | protected
        self.hide_emails = False

    @property
    def comments(self):
        if not self.comments:
            self.parse_comments()
        return self._comments

    @property
    def access(self):
        if not self.access:
            self.parse_comments()
        return self._access

    @property
    def arguments(self):
        raise NotImplementedError

    def parse_comments(self) -> None:
        """
        Fills the attributes ``.comments`` (dict), ``.arguments`` and ``access``.
        """
        # reset all temps
        self._previous_comments = ''
        self._ongoing_open_parenthesis = ''
        self._multiline = False
        self._current_class = None
        self._current_access = 'private'  # private | public | protected
        # parse
        filename = self.get_header_filename(self.target)
        with open(filename, 'r') as cpp_file:
            for line in cpp_file:
                line = line.strip()
                line = line.replace('enum class', 'class')
                self._parse_line(line)
        if self._previous_comments:
            self._comments['(EOF)'] += self._previous_comments
        # redact emails
        if self.hide_emails:
            for attribute, comment in self._comments.items():
                if re.search('[\w\._]+\@[\w\._]+\.\w+', comment):
                    self._comments[attribute] = re.sub('[\w\._]+\@[\w\._]+\.\w+', '👾👾.👾👾👾@👾👾👾.👾👾', comment)
        # done
        return self

    def _parse_line(self, line: str) -> None:
        if line == '':
            return
        elif self._multiline:
            if line.find('*/') == -1:
                self._previous_comments += line + '\n'
            else:
                self._multiline = False
                self._previous_comments += line.replace('*/', '') + '\n'
        # args are spanning multiple rows
        elif self._ongoing_open_parenthesis != '' and ')' not in line:
            self._ongoing_open_parenthesis += line
        elif self._ongoing_open_parenthesis != '':  # finished.
            self._ongoing_open_parenthesis += line.split(')')[0]
            self._parse_method_line(self._ongoing_open_parenthesis)
            self._ongoing_open_parenthesis = ''
        # crap comment lines
        elif line.find('// vi:') == 0:
            return
        elif line.find('/*') == 0:
            self._multiline = True
            self._previous_comments += line.replace('/*', '') + '\n'
        elif line.find('// (c)') == 0:
            return
        elif line.find('// -*-') == 0:
            return
        # good comment
        elif line.find('//') == 0:
            self._previous_comments += re.sub('/+', '', line) + '\n'
        # the top is finished.
        elif line.find('namespace') == 0:
            self._comments['__module__'] += self._previous_comments  # make it pythonic sounding...
            self._previous_comments = ''
        # its a class
        elif line.find('class ') == 0:
            # enum class will fail. hence earlier replace
            if re.match('class ([\w_]+)', line) is None:
                raise ValueError(f'line "{line}" failed parsing')
            self._current_class = re.match('class ([\w_]+)', line).group(1)
            self._current_access = 'private'  # default
            self._comments[self._current_class + '.__doc__'] += self._previous_comments
            self._previous_comments = ''
        elif line.find('public:') == 0:
            self._current_access = 'public'
        elif line.find('private:') == 0:
            self._current_access = 'private'
        elif line.find('protected:') == 0:
            self._current_access = 'protected'
        elif re.search('[\w_]+\s?\(.*\)', line):
            # something like void fold( core::pose::Pose& extended_pose, ProtocolOP prot_ptr );
            self._parse_method_line(line)
        elif re.search('[\w_]+\s?\(', line):
            self._ongoing_open_parenthesis += line

        else:
            return

    def _parse_method_line(self, line):
        # determine name of method
        parts = line.split('(')[0].strip().split()
        methodname = parts[-1]
        if methodname == self._current_class:
            method = '__init__'
        if self._current_class:
            fullname = self._current_class + '.' + methodname
        else:
            fullname = methodname
        # add comments
        if self._comments[fullname] != '':  # overloaded
            self._comments[fullname] += '.........'
        self._comments[fullname] += self._previous_comments
        # access
        for access in ('private', 'public', 'protected'):
            if access in parts:
                self._access[fullname] = access
                break
        else:
            self._access[fullname] = self._current_access
        self._previous_comments = ''

    def get_implemenatation_filename(self,
                                     target: Optional[pyrosetta.rosetta.protocols.moves.Mover] = None
                                     ) -> str:
        """
        C++ code is normally split into header files (.hh) or implementation files (.cpp/.cc).
        Namely, a class is declared in the former and the code is filled in the latter.
        """
        return self.get_src_filename(target, header=False)

    def get_header_filename(self,
                            target: Optional[pyrosetta.rosetta.protocols.moves.Mover] = None
                            ) -> str:
        """
        C++ code is normally split into header files (.hh) or implementation files (.cpp/.cc).
        Namely, a class is declared in the former and the code is filled in the latter.
        """
        return self.get_src_filename(target, header=True)

    def get_src_filename(self,
                         target: Optional[pyrosetta.rosetta.protocols.moves.Mover] = None,
                         header=True) -> str:
        if target is None:
            target = self.target
        if 'pyrosetta.rosetta.' not in target.__module__:
            raise ValueError(f'{target.__class__.__name__} is not a pybind11 Pyrosetta object')
        parts = target.__module__.replace('pyrosetta.rosetta.', '').split('.')
        mover_name = target.__class__.__name__
        if mover_name == 'pybind11_type':  # not an instance!
            mover_name = target.__name__
        if header is True:
            extensions = ('.h', '.hh')
        else:
            extensions = ('.cpp', '.cc', '.c')
        for extension in extensions:
            path = os.path.join(self.rosetta_folder, self.src_folder, *parts, mover_name + extension)
            if os.path.exists(path):
                return path
        else:
            raise FileNotFoundError(f'{path} does not exist')
