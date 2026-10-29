"""The simulation kernel as Python sees it.

Subdirectories follow the kernel's passes rather than its source files: the
junction, coding-end processing, D inversion, mutation and indels, productivity
assessment, receptor revision, and the corruption stack. Modules directly here
cover the engine surface itself -- pass plans, parameter signatures, recorded
provenance and replay.
"""
