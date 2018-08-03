def _init_table_h():
    _table_h = []
    for i in range(256):
        part_h = 0
        for j in range(8):
            rflag = i & 1
            i >>= 1
            if part_h & 1:
                i |= (1 << 31)
            part_h >>= 1
            if rflag:
                part_h ^= 0xd8000000
        _table_h.append(part_h)
    return _table_h


_TABLE_H = _init_table_h()


def crc64(s):
    """Returns the crc64 checksum for a sequence (string or Seq object).

    Copied from BioPython.

    Examples
    --------
    From UniParc:
    >>> crc64('MSGGKYVDSE')
    '368583B2DB533878'

    Note that the case is important:
    >>> crc64("ACGTACGTACGT")
    'C4FBB762C4A87EBD'
    >>> crc64("acgtACGTacgt")
    'DA4509DC64A87EBD'
    """
    if not s:
        return None
    _table_h = _TABLE_H  # speed up lookups by creating a local copy
    cdef long crcl = 0
    cdef long crch = 0
    cdef long shr, temp1h, temp1l, idx
    for c in s:
        shr = (crch & 0xFF) << 24
        temp1h = crch >> 8
        temp1l = (crcl >> 8) | shr
        idx = (crcl ^ ord(c)) & 0xFF
        crch = temp1h ^ _table_h[idx]
        crcl = temp1l
    return "%08X%08X" % (crch, crcl)
