"""Read selected columns of a remote Parquet file through HTTP range requests.

Only the footer and the requested column chunks are downloaded, so a single
column of a multi-GB public table can be scanned without fetching the whole file.
"""
import io
import os

import pyarrow.parquet as pq
import requests

CA = os.environ.get("REQUESTS_CA_BUNDLE") or os.environ.get("SSL_CERT_FILE") or True


class HTTPRangeFile(io.RawIOBase):
    def __init__(self, url, session=None, block=8 << 20):
        self.url, self.pos, self.block = url, 0, block
        self.s = session or requests.Session()
        r = self.s.head(url, timeout=60, verify=CA)
        r.raise_for_status()
        self.size = int(r.headers["Content-Length"])
        self._cache = {}
        self.bytes_fetched = 0

    def readable(self):
        return True

    def seekable(self):
        return True

    def tell(self):
        return self.pos

    def seek(self, offset, whence=0):
        self.pos = {0: offset, 1: self.pos + offset, 2: self.size + offset}[whence]
        return self.pos

    def _get(self, start, end):
        for attempt in range(4):
            try:
                r = self.s.get(self.url, headers={"Range": f"bytes={start}-{end - 1}"}, timeout=300, verify=CA)
                r.raise_for_status()
                self.bytes_fetched += len(r.content)
                return r.content
            except requests.RequestException:
                if attempt == 3:
                    raise

    def read(self, n=-1):
        if n is None or n < 0:
            n = self.size - self.pos
        n = min(n, self.size - self.pos)
        if n <= 0:
            return b""
        data = self._get(self.pos, self.pos + n)
        self.pos += len(data)
        return data

    def readinto(self, b):
        data = self.read(len(b))
        b[: len(data)] = data
        return len(data)


def read_columns(url, columns, filter_fn=None):
    """Return a pyarrow Table with `columns`, optionally filtered per row group."""
    f = HTTPRangeFile(url)
    pf = pq.ParquetFile(f)
    tables = []
    for rg in range(pf.num_row_groups):
        t = pf.read_row_group(rg, columns=columns)
        if filter_fn is not None:
            t = t.filter(filter_fn(t))
        tables.append(t)
    import pyarrow as pa
    return pa.concat_tables(tables) if tables else None, f.bytes_fetched
