"""
This module should be considered
temporary until all tools in here
are available via official release
channels.
"""

localrules: clone_nucfreqtwo_repo
rule clone_nucfreqtwo_repo:
    output:
        nucfreqtwo = "sources/NucFreqTwo/NucFreqTwo.py"
    envmodules:
        "git"
    params:
        repo_url = "git@git.hhu.de:ebertp/NucFreqTwo.git",
        branch_name = "split-two-phases"
    shell:
        "mkdir -p sources && cd source && "
        "git clone {params.repo_url} && "
        "cd NucFreqTwo && "
        "git switch {params.branch_name}"

