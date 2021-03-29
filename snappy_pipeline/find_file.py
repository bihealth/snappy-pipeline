# -*- coding: utf-8 -*-
"""Code for crawling the file system and caching the results
"""

import json
import logging
import os
import sys
from collections import OrderedDict
from fnmatch import fnmatch

from fasteners import InterProcessLock

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>"


class FileNamesTooDifferent(Exception):
    """Raised when two file names are too different to be PE reads"""


class PatternSet:
    """Store named or unnamed list of patterns"""

    def __init__(self, patterns, names=None):
        #: Patterns to search for with names
        self.patterns = tuple(patterns)
        #: Optional names
        self.names = tuple(names or [])
        if self.names and len(self.names) != len(self.patterns):
            raise ValueError(  # pragma: no cover
                "Names must be empty or have the same length as files ({} vs. {})".format(
                    self.names, self.patterns
                )
            )
        #: Named patterns, if any, else ``None``
        self.named_patterns = None
        if self.names:
            self.named_patterns = OrderedDict(zip(names, patterns))

    def __str__(self):
        return "PatternSet({}, {})".format(self.patterns, self.names)

    def __repr__(self):
        return str(self)


class FileSystemCrawlerResult:
    """n-tuple of optionally named files"""

    def __init__(self, base_folder, files, names=None):
        #: Folder to start crawling in
        self.base_folder = base_folder
        #: Patterns to search for
        self.files = tuple(files)
        #: Names for the file patterns, optional; if given has to have the same length as files
        self.names = tuple(names or [])
        if self.names and len(self.names) != len(self.files):
            raise ValueError(  # pragma: no cover
                "Names must be empty or have the same length as files ({} vs. {})".format(
                    self.names, self.files
                )
            )
        #: Dict with name-to-pattern mapping, ``None`` if ``names`` is not given
        self.named_files = None
        if self.names:
            self.named_files = OrderedDict(zip(names, files))

    def to_dict(self):
        """Convert to dict, can only work if self.names and self.files is given"""
        # TODO: remove?
        if not self.names:
            raise ValueError("No names, cannot convert to dict")
        return dict(zip(self.names, self.files))

    def __str__(self):
        tpl = "FileSystemCrawlerResult({})"
        return tpl.format(
            ", ".join(map(repr, (self.base_folder, self.files, self.names, self.named_files)))
        )

    def __repr__(self):
        return str(self)


class FileSystemCrawler:
    """Crawl the file system

    - start crawling the file system from a given directory
    - look for files matching a given ``PatternSet``
    - that are below a directory with a given name
    """

    cache_version = 1

    def __init__(self, cache_path, invalidation_paths, lock_timeout=60):
        #: The logger to use.
        self.logger = logging.getLogger("file_crawler")
        #: Path to cache (will be stored in JSON format)
        self.cache_path = cache_path
        #: Path to files to use for checking invalidation.
        self.invalidation_paths = invalidation_paths
        #: The actual dict with the cache, loaded from path to ``cache_path`` if the cache file
        #: exists.
        self.cache = None
        #: Flag whether cache has been modified and needs saving
        self.cache_dirty = False
        #: Flag whether cache has been invalidated already.
        self.cache_invalidated = False
        #: Timeout for obtaining file system lock on the file system
        self.lock_timeout = lock_timeout
        if os.path.exists(self.cache_path):
            self.cache_invalidated = False
            self.cache = self._load_cache()
        else:
            self._set_fresh_cache()

    def _set_fresh_cache(self):
        """Set cache to a fresh state."""
        self.cache_invalidated = True
        self.cache_dirty = True
        self.cache = {"cache_version": self.__class__.cache_version, "root_dirs": {}}

    def run(self, root_dir, dir_name, pattern_sets, allow_empty_right):
        """Perform the file system crawling from a root directory given a query pattern set

        ``allow_empty_right`` -- for mixed PE/SE read data sets (must be either SE or PE
                                 for one library!)
        """
        matches = {}  # path => pattern set idx => pattern idx => [path]
        returned = 0  # from {0, 1, 2}; how many patterns matched?
        # Invalidate cache fully if the cache file is older than any one of self.invalidation_paths
        self._perform_cache_invalidation()
        # Ensure that cache entry with crawling results of all files exists
        if root_dir not in self.cache["root_dirs"]:
            self.cache_dirty = True
            self.cache["root_dirs"][root_dir] = tuple(sorted(self._build_cache(root_dir)))
        # Now, crawl over this structure and match against all pattern sets
        self.logger.debug('Crawling "%s" for dir_name "%s"', root_dir, dir_name)
        for i, pattern_set in enumerate(pattern_sets):
            self.logger.debug("  patterns in pattern set #%d: %s", i, pattern_set.patterns)
        for path in self.cache["root_dirs"][root_dir]:
            if dir_name not in path:
                continue
            idx = path.index(dir_name)
            left, right = path[: idx + 1], path[idx + 1 :]
            for i, pattern_set in enumerate(pattern_sets):
                for j, pattern in enumerate(pattern_set.patterns):
                    does_match = fnmatch("/".join(right), pattern)
                    self.logger.debug(
                        'does "%s" match "%s" match? => %s', "/".join(right), pattern, does_match
                    )
                    if does_match:
                        matches.setdefault("/".join(left), {}).setdefault(i, {}).setdefault(
                            j, []
                        ).append("/".join(right))
        # Go over results and check whether they are conforming.
        for path, path_matches in matches.items():
            for set_idx, set_matches in path_matches.items():
                # Must have a match for each pattern
                if not allow_empty_right and len(set_matches) != len(
                    pattern_sets[set_idx].patterns
                ):
                    print(
                        (
                            "WARNING: Skipping matches {} as the number of matches is not equal to "
                            "the number of patterns in {}"
                        ).format(set_matches, pattern_sets[set_idx].patterns),
                        file=sys.stderr,
                    )
                    continue
                # Guard against mixing SE and PE results for crawling
                if returned:
                    if returned != len(set_matches):
                        raise ValueError(  # pragma: no cover
                            "Found mixed SE and PE data for one library!"
                        )
                else:
                    returned = len(set_matches)
                # Must have the same number of matches for each pattern
                lst_lens = [len(l) for l in set_matches.values()]
                if len(set(lst_lens)) != 1:
                    raise ValueError(  # pragma: no cover
                        "Must have the same number of matches per pattern, but found {}".format(
                            list(set_matches.values())
                        )
                    )
                # Yield result, checking that file names are equal except for one character (e.g.,
                # "R1" vs "R2")"
                for i in range(0, lst_lens[0]):
                    files = [
                        os.path.join(root_dir, path, match[i]) for match in set_matches.values()
                    ]
                    self._check_fname_mismatches(files)
                    base_path = os.path.join(root_dir, path)
                    yield FileSystemCrawlerResult(
                        base_path, files, pattern_sets[set_idx].names[:returned]
                    )

    def _perform_cache_invalidation(self):
        """Check whether the cache needs to be invalidated and do so if necessary."""
        if self.cache_invalidated:
            return  # Cache has been invalidated before
        if not os.path.exists(self.cache_path):
            return  # No cache yet
        self.logger.debug("Checking whether file crawler cache should be invalidated...")
        cache_ctime = os.path.getctime(self.cache_path)
        for path in self.invalidation_paths:
            path_mtime = os.path.getmtime(path)
            if path_mtime > cache_ctime:
                self.logger.info("Invalidating cache because of %s", path)
                self._set_fresh_cache()
                return
        self.logger.debug(" => no, not invalidating cache")

    @classmethod
    def _check_fname_mismatches(cls, file_names):
        for i in range(0, len(file_names)):
            for j in range(0, i):
                a = file_names[i]
                b = file_names[j]
                mismatches = 0
                if len(a) != len(b):
                    raise FileNamesTooDifferent(  # pragma: no cover
                        "File names have different length {} vs {}".format(a, b)
                    )
                for x, y in zip(a, b):
                    mismatches += int(x != y)
                if mismatches > 1:
                    raise FileNamesTooDifferent(  # pragma: no cover
                        "File names too different ({} mismatches) {} vs {}".format(mismatches, a, b)
                    )

    def _build_cache(self, root_dir):
        self.logger.info("Building file system crawler cache from %s", root_dir)
        for root, _, files in os.walk(root_dir, followlinks=True):
            self.logger.debug("Caching for directory %s", root)
            base = root[len(root_dir) + 1 :].split("/") or ["."]
            yield from (tuple(base + [f]) for f in files)

    def save_cache(self, cache_path=None):
        """Save cache, ``cache_path`` overriding ``self.cache_path``"""
        if not (self.cache_dirty or self.cache_invalidated):
            return  # don't save if unchanged
        cache_path = cache_path or self.cache_path
        with InterProcessLock(self.cache_path + ".lock"):
            self.logger.debug("Saving file system crawler cache to %s", cache_path)
            with open(cache_path, "wt") as f:
                json.dump(self.cache, f)

    def _load_cache(self):
        with InterProcessLock(self.cache_path + ".lock"):
            self.logger.info("Loading file system crawler cache from %s", self.cache_path)
            with open(self.cache_path, "rt") as f:
                result = json.loads(f.read(), object_pairs_hook=OrderedDict)
        if result["cache_version"] != self.__class__.cache_version:
            raise ValueError(  # pragma: no cover
                "Invalid cache version {}".format(result["cache_version"])
            )
        return result
