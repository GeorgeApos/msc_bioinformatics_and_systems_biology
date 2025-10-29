#!/usr/bin/env python3
"""
Author: Sandra Smit
Script to print the user's name and the python version
"""

from sys import version
from os import getcwd


def main():
    """Main function of this script
    """

    name = ""

    # This line generates output, don't touch it
    print("Hello {}!".format(name))
    print(getcwd())
    print("You are running this script with Python version: ", version)


if __name__ == "__main__":
    main()