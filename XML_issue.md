## Issue encountered

One problem I faced with the XML parser is that I am not too sure how to use `-parser:script_vars` successfully.

This is compounded by the fact I am not sure how to use the `RosettaScriptsParser` properly in the first place.

Which in turn is compounded by the fact I have no clue about `OptionCollection`...

These are my notes so future-me can reference them, someone else may learn something or someone may know the solution...

## RosettaScriptsParser

This object has lots of fancy methods to generate a protocol:

    xml = pyrosetta.rosetta.protocols.rosetta_scripts.RosettaScriptsParser()

A protocol is a mover that runs like a regular mover, `protocol.apply(pose)` and is basically the PROTOCOL tag of the XML.
Below `generate_mover_and_apply_to_pose` gets called in the tests. 
This is accepts a bare minimum of arguments —no blasted `OptionCollection` as discussed below.
It says "apply_to_pose" but I don't thing it finished by calling `protocol.apply(pose)`.

A lot of operations it does involve tags, `pyrosetta.rosetta.utility.tag.Tag` objects:

    tag = xml.create_tag_from_xml(xml_fname=xml_filename)

This does not do much in itself, but can be passed to the xml parser:

    xml.generate_mover_for_protocol(pose=pyrosetta.Pose(),
                                    modified_pose=False,
                                    protocol_tag=tag,
                                    options=options)
                                    
If you can get the `OptionCollection` object `options` to work...

## OptionCollection and -parser:script_vars

Several objects require a `OptionCollection`.
This can be instantiated normally:

    options = pyrosetta.rosetta.utility.options.OptionCollection()
    
But this when used:

    xml = pyrosetta.rosetta.protocols.rosetta_scripts.RosettaScriptsParser()
    tag = xml.create_tag_from_xml(xml_fname='/well/brc/matteo/rosetta_bin_linux_2020.08.61146_bundle/main/tests/integration/tests/farnesyl/farnesyl.xml',
                           script_vars=pyrosetta.rosetta.utility.vector1_std_string())
    options = pyrosetta.rosetta.utility.options.OptionCollection()
        
    xml.generate_mover_for_protocol(pose=pyrosetta.Pose(),
                                    modified_pose=False,
                                    protocol_tag=tag,
                                    options=options)
                                
results in the error you get when there is a command line key missing:

    File: /scratch/benchmark/W.hojo-2/rosetta.Hojo-2/_commits_/main/source/build/PyRosetta/linux/clang-3.4.2/python-3.7/release/source/src/utility/keys/SmallKeyVector.hh:538
    [ ERROR ] UtilityExitException
    ERROR: Assertion `active( key )` failed.

Some functions return this `OptionCollection`...

    options = pyrosetta.rosetta.basic.options.initialize()
    
However, this will result in a segfault.
Stuff can be added to an `OptionCollection` initialised like above:

    options.load_options_from_file(flags_filename)

While calling the class constructor initilialised version will result in `-in:file:s` is unknown...

What is peculiar is that I cannot seem to get the values off it.

    ok = pyrosetta.rosetta.utility.options.StringVectorOptionKey('-parser:script_vars')
    options.has(ok) # False
    
Yet `print(options)` lists it.

There are a few `show` methods, but only `show_user` gives the non-generic results, but the same as print.

    buffer = pyrosetta.rosetta.std.stringbuf()
    options.show_user(pyrosetta.rosetta.std.ostream(buffer))
    print(buffer.str().strip())

However, the sole reason to parse it is to get `script_vars`, which I can extract myself:

    import re
    flags_filename = '...'
    flags = pyrosetta.rosetta.utility.vector1_std_string()
    flags.extend( re.findall('-parser:script_vars ([^\n\-]+)', open(flags_filename).read()) )

However, whereas this may work on a xml script without script_vars:

    xml = pyrosetta.rosetta.protocols.rosetta_scripts.RosettaScriptsParser()
    xml.generate_mover_and_apply_to_pose(pose=pyrosetta.Pose(),
                                         xml_fname=xml_filename)

The version with script_vars will segfault:

    xml = pyrosetta.rosetta.protocols.rosetta_scripts.RosettaScriptsParser()
    xml.generate_mover_and_apply_to_pose(pose=pyrosetta.Pose(),
                                         modified_pose=False,
                                         script_vars=flags,
                                         xml_fname=xml_filename)
                                         
Damn. Maybe the method called is wrong?

