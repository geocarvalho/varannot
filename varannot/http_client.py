"""
http_client.py
==============
Shared HTTP client with on-disk caching and polite rate limiting.

All API responses are cached as JSON files keyed by a hash of the request,
so re-running the tool on the same variants does not re-hit the network.
"""

import hashlib
import json
import os
import time

import requests


class CachedSession:
    """A requests session wrapper with file-based caching and rate limiting."""

    def __init__(self, cache_dir=".varannot_cache", min_interval=0.4, timeout=30):
        self.cache_dir = cache_dir
        self.min_interval = min_interval   # seconds between live requests
        self.timeout = timeout
        self._last_request = 0.0
        os.makedirs(cache_dir, exist_ok=True)
        self.session = requests.Session()
        self.session.headers.update({"User-Agent": "varannot/1.0 (research use)"})

    def _cache_path(self, key):
        h = hashlib.sha256(key.encode("utf-8")).hexdigest()[:32]
        return os.path.join(self.cache_dir, f"{h}.json")

    def _throttle(self):
        elapsed = time.time() - self._last_request
        if elapsed < self.min_interval:
            time.sleep(self.min_interval - elapsed)
        self._last_request = time.time()

    def get_json(self, url, params=None, headers=None, cache_key=None, timeout=None):
        """GET request returning parsed JSON, with caching."""
        key = cache_key or f"GET:{url}:{json.dumps(params, sort_keys=True)}"
        path = self._cache_path(key)
        if os.path.exists(path):
            with open(path, "r", encoding="utf-8") as fh:
                return json.load(fh)

        self._throttle()
        merged_headers = {"Accept": "application/json"}
        if headers:
            merged_headers.update(headers)
        cacheable = True
        try:
            resp = self.session.get(url, params=params, headers=merged_headers,
                                    timeout=timeout or self.timeout)
            if resp.status_code == 200:
                data = resp.json()
            elif resp.status_code in (400, 404):
                data = {"_error": f"HTTP {resp.status_code}", "_status": resp.status_code}
            else:
                resp.raise_for_status()
                data = resp.json()
        except (requests.RequestException, ValueError) as exc:
            # Network/timeout/parse failures are transient: don't cache them,
            # so a later run can retry instead of being stuck on the error.
            data = {"_error": str(exc)}
            cacheable = False

        if cacheable:
            with open(path, "w", encoding="utf-8") as fh:
                json.dump(data, fh)
        return data

    def post_json(self, url, payload, headers=None, cache_key=None, timeout=None):
        """POST request returning parsed JSON, with caching."""
        key = cache_key or f"POST:{url}:{json.dumps(payload, sort_keys=True)}"
        path = self._cache_path(key)
        if os.path.exists(path):
            with open(path, "r", encoding="utf-8") as fh:
                return json.load(fh)

        self._throttle()
        merged_headers = {"Accept": "application/json", "Content-Type": "application/json"}
        if headers:
            merged_headers.update(headers)
        cacheable = True
        try:
            resp = self.session.post(url, json=payload, headers=merged_headers,
                                     timeout=timeout or self.timeout)
            if resp.status_code == 200:
                data = resp.json()
            elif resp.status_code in (400, 404):
                data = {"_error": f"HTTP {resp.status_code}", "_status": resp.status_code}
            else:
                resp.raise_for_status()
                data = resp.json()
        except (requests.RequestException, ValueError) as exc:
            data = {"_error": str(exc)}
            cacheable = False

        if cacheable:
            with open(path, "w", encoding="utf-8") as fh:
                json.dump(data, fh)
        return data
