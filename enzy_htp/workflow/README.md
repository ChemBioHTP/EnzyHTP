# Workflow

Author: Zhong, Yinjie.  
Email: yinjie.zhong@vanderbilt.edu

[TOC]

## 1. Introduction

This module is to assemble JSON format task configuration files written in the specified format into a computing task, and supports saving the task as a binary file, re-reading the binary file, and reloading the task by reading the updated configuration file after the task exits with an error, as well as to continue execution from the updated/error-reporting position.

## 2. Quick Start

In this part, we will be provided with a step-by-step manual to learn the basic usage of this module, so as to be able to run our first computation job in 10 minutes.

### 2.1 Simple example

Before initializing your job, we need to prepare a json format configuration file as an _order_ (just like making an order of dishes you want to a restaurant).

Here's an example of a JSON configuration file (hereafter referred to as an **order**)

```json
{
  "api": "general",
  "store_as": "general",
  "args": {
    "workunits": [
      {
        "api": "read_pdb",
        "store_as": "structure",
        "args": {
          "path": "./data/7si9_rm_water_disconnected.pdb"
        }
      },
      {
        "api": "remove_hydrogens",
        "store_as": "structure",
        "args": {
          "stru": "`structure`",
          "polypeptide_only": true
        }
      },
      {
        "api": "protonate_stru",
        "args": {
          "stru": "`structure`",
          "ph": 7.0,
          "protonate_ligand": false,
          "engine": "pdb2pqr"
        }
      }
    ]
  }
}
```

In an order, each pair of large bracket (`{` and `}`) containing keys like `api`, `store_as` and `args` represents a workunit, and each pair of middle bracket (`[` and `]`) containing several workunits represents a workflow.

In each workunit,

- The value of the `api` key indicates the function called by the workunit and determines the type and function of the workunit;
- The value of the `store_as` key declares the name of the variable into which the return value of the function called will be stored;
- The value of the `args` key stores the values of the arguments that are to be used for the called function in the form of a dictionary.

In this order, the outermost level is a workunit with an `api` of `general`, which is a workunit for process control purposes, and it is only by using it as the outermost level that we can use the full features of this module.

The other `api` values we see, such as `read_pdb`, `remove_hydrogens`, `protonate_stru`, etc., are specific Science API values that are mapped to a callable object (function).

The way to read and run this order is as follows.

```python
from enzy_htp.workflow import ExecutionEntity, WorkFlow, WorkUnit, GeneralWorkUnit
from enzy_htp.workflow import SCIENCE_API_MAPPER, StatusCode

general = GeneralWorkUnit.from_json_filepath(json_filepath="path/to/order.json", working_directory="your/working/directory", save_snapshot=True, overwrite_database=True)
return_key, return_value = general.execute()
pickle_filepath = general.latest_pickle_filepath
```

When you created the `general` in the code above, a SQLite database file is created in our working directory (as you specified, using our current directory as default), which will update the running status of your.

The data output throughout the execution of the `general` are able to be read from the `return_value`.

A snapshot is saved as a binary file to the `pickle_filepath` file at the end of execution.

- If the order we submit runs fine the first time we run it, we don't need to make any updates and just go ahead and read the data it outputs; 
- If we face errors while running it, and you need to update your order but don't want to run it from scratch (because some of the steps are very time-consuming), you can use the following method to do so.

To commit any updates to your task, just make modifications to your order, and run following commands.

```python
reloaded_general = GeneralWorkUnit.load_snapshot_file(filepath=pickle_filepath)
reloaded_general.reload(json_filepath="path/to/updated_order.json")
reloaded_general.execute()
```

Here, I'm using the `reloaded_general` variable to make it easier to distinguish between reloaded `general`, whereas you could actually just use the `general` variable and have the old object overwritten after executing `GeneralWorkUnit.load_snapshot_file`.

### 2.2 Assign values

Assign values at the outermost level (The `general` level) is more recommended due to convenience...

(To be continued...)

### 2.3 Loop Structure

Descriptions about loop structure...

(To be continued...)

## 3. Advanced Manual