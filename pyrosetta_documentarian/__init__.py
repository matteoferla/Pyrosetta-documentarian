from .attributes import AttributeDocumentarian  # document what the values of the attributes are
from .comment import CommentDocumentarian  # get the C++ comments
from .xml import XMLDocumentarian  # find the relevant XML script test

class Documentarian(AttributeDocumentarian, XMLDocumentarian, CommentDocumentarian):
    # ``CommentDocumentarian`` has its own init.
    pass