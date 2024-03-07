#!/usr/bin/env python
import sys

operators = {
    'move':       0b000_00000,
    'mid_move':   0b100_00000_00000000,
    'long_move':  0b111_00000_00000000_00000000_00000000,
    'change':     0b001_00000,
    'mid_change': 0b011_00000_00000000,
    'next':       0b010_00000_00000000_00000000_00000000,
    'zero':       0b101_00000,
    'chmo':       0b111_00000,
}

MAX_SHORT = 2**5
MAX_MID = 2**13
MAX_LONG = 2**29


class Code:
    def __init__(self):
        self.byte_list = []
        self.last_index_point = 0
        self.index_points = []

    def add(self, code_point, contig=None, pos=None):
        filepos = len(self.byte_list)
        if code_point.kind == "next" or (code_point.kind == "change" and filepos - self.last_index_point.filepos > 10000):
            ip = IndexPoint(filepos, contig, pos)
            self.index_points.append(ip)
            self.last_index_point = ip
        self.byte_list.extend(code_point.byte_list())


    def next_contig(self, dp, contig=None, pos=None):
        assert dp <= MAX_LONG, f"Too high initial DP for new contig: {dp}"
        self.add(CodePoint("next", dp), contig, pos)

    def zero(self, n):
        self.add(CodePoint("zero", min(n, MAX_SHORT)))
        if n > MAX_SHORT:
            self.move(n-32)

    def chmo(self, n, bp):
        assert bp <= 8 and n + 2 <= 4, f"Cannot make change-move with n={n}, bp={bp}"
        if n < 0:
            value = (n + 2) + bp * 4
        else:
            value = (n + 1) + bp * 4

        self.add(CodePoint("chmo", n))

    def change(self, n, contig, pos):
        assert n <= MAX_MID, f"Too large change of DP: {n}"
        if n+16 <= MAX_SHORT and n+16 >= 0:
            self.add(CodePoint("change", n+16), contig=contig, pos=pos)
        else:
            self.add(CodePoint("mid_change", n + MAX_MID//2), contig=contig, pos=pos)

    def move(self, bp):
        assert bp <= MAX_LONG, f"Too long move in bp: {bp}"
        if bp <= MAX_SHORT:
            self.add(CodePoint("move", bp))
        elif bp <= MAX_MID:
            self.add(CodePoint("mid_move", bp))
        else:
            self.add(CodePoint("long_move", bp))


    def write_to_file(self, filename):
        with open(filename, "wb") as out_fh:
            out_fh.write(bytearray(self.byte_list))

class IndexPoint:
    def __init__(self, filepos, contig, pos):
        self.filepos = filepos
        self.contig = contig
        self.pos = pos

    def __repr__(self):
        return f"{self.filepos}\t{self.contig}:{self.pos}"

class CodePoint:
    def __init__(self, kind, value):
        self.kind = kind
        self.value = value
        if kind in ["move","change","zero", "chmo"]:
            self.nbytes = 1
        elif kind in ["mid_move", "mid_change"]:
            self.nbytes = 2
        elif kind in ["long_move", "next"]:
            self.nbytes = 4
        else:
            print("Unknown codepoint kind")
            sys.exit(1)


    def binary(self):
        return bin(operators[self.kind] + self.value)

    def byte_list(self):
        bl = []
        if self.nbytes == 1:
            bl.append(operators[self.kind] + self.value)
        if self.nbytes == 2:
            i16 = operators[self.kind] + self.value
            bl.append( (i16 >> 8) & 0xFF )
            bl.append(i16 & 0xFF)
        if self.nbytes == 4:
            i32 = operators[self.kind] + self.value
            bl.append( (i32 >> 24) & 0xFF )
            bl.append( (i32 >> 16) & 0xFF )
            bl.append( (i32 >> 8) & 0xFF )
            bl.append(i32 & 0xFF)

        return bl


    def __repr__(self):
        val = self.value
        if self.kind == "change":
            val -= 16
        string = f"{self.kind}\t{val}\t{self.nbytes}\t{self.binary()}\t"
        for b in self.byte_list():
            string += str(b)+","
        return string


def main():
    code = Code()
    with open(sys.argv[1]) as fh:
        prev_contig, prev_dp, prev_end = None, 0, 0
        for idx, line in enumerate(fh):
            contig, start, end, dp = line.strip().split("\t")
            end, start, dp = int(end), int(start), int(dp)
            region_len = end - start

            if not prev_contig or contig != prev_contig:
                code.next_contig(dp, contig, start)
                code.move(region_len)
            elif start != prev_end:
                code.zero(start-prev_end)
            else :
                dp_change = dp - prev_dp
                if abs(dp_change) <= 2 and region_len <= 8:
                    code.chmo(dp_change, region_len)
                else:
                    if dp_change != 0:
                        code.change(dp_change, contig, start)
                    code.move(region_len)

            prev_dp, prev_contig, prev_end = dp, contig, end


    for idx, index_point in enumerate(code.index_points):
        print(idx, index_point)

    code.write_to_file("out.bd")

if __name__ == "__main__":
    main()
