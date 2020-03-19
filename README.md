# AGR File Generator

## Description

This tool creates files from Alliance resources. The tool is build and testsed with python 3

## Current File Types

- VCF
- DAF
- Orthology
- Expression

## Running Tool Example command

You will want to customize your run by setting these two env variable in the docker-compose.yaml file

```bash
 export ALLIANCE_RELEASE='2.3.0'
 export NEO4J_HOST=build.alliancegenome.org
```

if you only want to generate certain types of files make sure to comment/uncomment the files in the src/agr/app.py file.

## Generate container

```bash
 make build
```

## generate files

```bash
 make run
```

## Generated Files

Files that are generated are located in the docker shared folder. This folder can be found by using the following command:

```bash
 docker inspect <container>
```

## Run tests

```bash
 make test
```

##Venv

```bash
source venv/bin/activate
```
