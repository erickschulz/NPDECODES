#! /usr/bin/env python
# -*- coding: utf-8 -*-

import os
import json
import subprocess
import tempfile
import zipfile
import shutil


def copy(file, indir, outdir):
    """
    Aggressively copy file to outdir. File can be a directory and can
    already exists.

    :param file: file (or directory) name
    :param indir: the directory which contains file
    :param outdir: the destination to where file copies to
    :return:
    """
    print("Moving '{}' to '{}'".format(indir + file, outdir + file))
    if file == "":
        shutil.copytree(indir, outdir)
        # Recursively copy an entire directory tree rooted at indir
        # to a directory named outdir and return the destination directory.
    elif file[-1] == "/":
        # if the solution needed to copy is a directory, delete the same directory in outdir (if has)
        # then use copytree to copy the directory
        try:
            shutil.rmtree(outdir + file)  # Delete an entire directory tree
        except:
            pass
        shutil.copytree(indir + file, outdir + file)
    else:
        # copy the file
        shutil.copy2(indir + file, outdir + file)
        """
        shutil.copy(src, dst, *, follow_symlinks=True)
        Copies the file scr to the file or directory dst. src and dst should be strings. If dst specifies a directory, 
        the file will be copied into dst using the base filename from src. Returns the path to the newly created file.
        """


def deploy(filename, indir, outdir,
           with_solution=False, with_internal=False):
    """
    Parse a C++ file trough "unifdef" and remove SOLUTIONS and INTERNAL ifdef
    code blocks. Must have "unifdef" program. Can toggle if INTERNAL and SOLUTIONS are
    true or false.

    :param filename:  e.g. "1dwaveabsorbingbc.cc"
    :param indir:     "../homeworks/1DWaveAbsorbingBC/mastersolution"
    :param outdir:    "../homeworks/1DWaveAbsorbingBC/templates"
    :param with_solution: remove block in #if SOLUTION
    :param with_internal:
    :return:
    """
    if is_cpp(indir + filename) or is_hpp(indir + filename):
        cmd = ["unifdef"]
        if with_internal:
            cmd.append("-DINTERNAL=1")
        else:
            cmd.append("-DINTERNAL=0")
        if with_solution:
            cmd.append("-DSOLUTION=1")
        else:
            cmd.append("-DSOLUTION=0")

        cmd.append(indir + filename)

        print("Preparing templates for '{}'".format(filename))
        print(" ".join(cmd))

        f = open(outdir + filename, "w")
        subprocess.call(cmd, stdout=f)  # run the command described by cmd
        f.close()
    else:
        copy(filename, indir, outdir)


def mkdir(dir):
    """
    Aggressively create a directory if non existing.
    """
    if os.path.exists(dir):
        shutil.rmtree(dir)
    os.makedirs(dir)


def is_cpp(file):
    """
    Returns True if file looks like a C++ file (header of .cpp)
    """
    # return "cpp" == file.split(".")[-1]
    return file.split(".")[-1] in ["cpp", "cc", "c"]


def is_hpp(file):
    """
    Returns True if file looks like a C++ file (header of .cpp)
    """
    return file.split(".")[-1] in ["hpp", "h"]


def generate_templates(problem_dir):
    """
    First create templates directory
    then copy the cpp/hpp files in mastersolution to templates folder

    :param problem_dir: "../homeworks/problem_name/"
    """
    mkdir(problem_dir + "templates/")  # "../homeworks/problem_name/templates"
    fileInMastersolution = os.listdir(problem_dir + 'mastersolution/')
    for file in fileInMastersolution:
        file_path = problem_dir + 'mastersolution/' + file
        if is_cpp(file_path) or is_hpp(file_path):
            deploy(file, problem_dir + 'mastersolution/', problem_dir + 'templates/', False)


def parse_json(filename):
    """
    Parse a JSON description of the Assignment bundles.
    """
    this_dir = os.path.dirname(os.path.realpath(__file__))
    os.chdir(this_dir)  # change the current working directory to this_dir

    f = open(filename, 'r')

    obj = json.load(f)

    assignment_dir = obj["assignment_dir"]  # "../Assignments/Codes"
    working_dir = obj["working_dir"]  # "./" current directory
    if not assignment_dir.endswith("/"):
        assignment_dir += "/"
    if not working_dir.endswith("/"):
        working_dir += "/"

    print("Looking for files in '{}'".format(assignment_dir))

    for problem in obj["ProblemSheets"]:
        problem_dir = assignment_dir + problem
        generate_templates(problem_dir)


if __name__ == "__main__":
    try:
        subprocess.call(["unifdef", "--help"])
    except FileNotFoundError:
        print("You must install 'unifdef'.")
        exit(-1)

    parse_json("assignment_list.json")

