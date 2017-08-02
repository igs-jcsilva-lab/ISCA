# Running the tool via Common Workflow Language (CWL)

CWL Documentation can be found [here](http://www.commonwl.org/draft-3/UserGuide.html#Writing_Workflows).

## Setup
1. Modify `conf.py` to have the correct paths for the installations of the various dependencies
2. Run the `setup.py` script to modify the CWL files to correctly utilize these installations
3. Setup complete, can move on to [Invoking a CWL workflow](https://github.com/jmatsumura/targeted_assembly/tree/master/cwl#invoking-a-cwl-workflow)

TIP: If you mess up a path in the conf and need to remake your __.cwl__ files, backup files are made by default (with extension of __.bak__) to allow one to rerun the `setup.py` script. Use the following command once you are in the same directory as all of your __.bak__ files, then go back to step 1 of [Setup](https://github.com/jmatsumura/targeted_assembly/tree/master/cwl#setup).
```
find . -iname "*.bak" -exec bash -c 'mv "$0" "${0%\.cwl\.bak}.cwl"' {} \;
```

## Invoking a CWL workflow
```
cwl-runner <cwl tool/workflow script> <input parameter yml/json>
```
The first parameter is a valid cwl tool or workflow script.  These have the extension __.cwl__.

The second parameter is a YAML or JSON file consisting of input parameters for the CWL script. YAML examples are provided and are listed with the extension __.yml__.

For example, to run the complete workflow/pipeline __first modify__ `targeted_asembly.yml`. Be sure to read and change as appropriate any line that does __not__ start with a __#__. Also make sure to have the [required inputs](https://github.com/jmatsumura/targeted_assembly#required-inputs) ready as they must be set in this file. Once `targeted_assembly.yml` has been modified appropriately, use the command:
```
cwl-runner --outdir /path/to/my/results targeted_assembly.cwl targeted_assembly.yml
```
`targeted_asembly.cwl` is comprised of all the other __.cwl__ files in this directory and will run them in the correct order. This command will output all the results of the workflow to the directory specified after the `--outdir` option. Output defaults to the current working directory if the `--outdir` option is not specified. 

### Outputs
- The primary output file will always be named `assembled_seqs.fsa` and consists of all the assembled sequences in a single file. 
- It is possible that all sequences will be assembled before the workflow reaches its end. In some cases, it may appear that the workflow has failed when all that has happened is there is nothing left to assemble. At each relevant step of the workflow, output files will be generated at the declared `--outdir`. Thus, in the event of some failures (those not related to erroneous inputs/syntax), one can check for the presence of any `*sequences.fsa` files that were successfully attained along the way and see if it comprises their complete data set. 

## Testing/Debugging
Example __.yml__  files which correspond to each __.cwl__ file can be found in the `individual_component_yml_examples` subdirectory. Modification of these to fit one's data allows one to run each step of the pipeline on its own for debugging purposes. 
