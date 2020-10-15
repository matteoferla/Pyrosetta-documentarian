# Pyrosetta-documentarian
A class to help reverse engineer what what a Pyrosetta class does...

## Note

I have combined various codeblocks I have written in the past to figure out how a Pyrosetta class works.
I plan to run it through all the classes and compile a set of tables of results.

## Usage

Import it

    from pyrosetta_documentarian import Documentarian
    Documentarian.rosetta_folder = 'User/matteo/rosetta'
    
That is, to use the XML and the Comment functionalities it needs a regular Rosetta path to work.

Given a mover

    mover = ...
    
get information for it

    doco = Documentarian(mover) # doco.target == mover
    print(doco.citation) # does `doco.target` have a citation on record?
   
### Find relevant XML test script
 
Here is the functionality from the `XMLDocumentarian` base class:
    
    # Get the paths of the XML test scripts that use the given mover.
    paths = doco.get_relevant_scripts() 
    # this takes a while as it fills .mover_directory by reading all the scripts.
    
    assert paths, 'No script uses it'
    reference = doco.get_mover_from_script(paths[0])

### See differing attributes
    
Here is the unctionality from the AttributeDocumentarian base class
    
    doco.compare(reference)

The latter returns a pandas table showing what the differences are between `doco.target` and `reference`.
Alternatively,

    doco.compare()

Would compare `doco.target` to a blank instance of its parent class (`doco.target_cls`).
    
### See what are the C++ comments
The documentation of pyrosetta is derived from `@brief` tags in comments preceeding methods and classes.
This assumes there is nothing odd.
Additionally due to an annoying feature of Sphynx the documentation skips `__init__` where most info is stored.

Here is the functionality from the CommentDocumentarian base class

    doco.hide_emails = True
    print(doco.comments)

