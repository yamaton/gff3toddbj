#!/usr/bin/env python
"""
Unquote unicode characters


Example:
From: Note=[KO:K20221]%3Babbreviation%3DIPO4%2C RANBP4%3Bsynonyms%3Dimportin-4
To  : Note=[KO:K20221];abbreviation=IPO4, RANBP4;synonyms=importin-4

"""

import fileinput
import urllib.parse

if __name__ == "__main__":
    with fileinput.input(mode="r") as f:
        for line in f:
            print(urllib.parse.unquote(line), end='')
