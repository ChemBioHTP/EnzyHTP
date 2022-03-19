from typing import Any


class MMPBSAConfig:
    def __init__(self, parent=None):
        self._parent = parent

    # -----------------------------
    # User defined MMPBSA.py(.MPI) exe dir other than
    #
    HOME = "$AMBERHOME"
    EXE = "$AMBERHOME/bin/MMPBSA.py.MPI"

    def required_executables(self):
        return [self.EXE]

    def required_env_vars(self):
        return [self.HOME]

    # ==============================
    # Method to express MMPBSA EXE path
    #
    @classmethod
    def get_MMPBSA_engine(cls):
        """
        Give default value to MMPBSA engine
        Only support MPI version now.
        ---
        return engine_path
        """
        if cls.MMPBSA_EXE == None:
            engine_path = Config.Amber.AmberHome + "/bin/MMPBSA.py.MPI"
        else:
            engine_path = cls.MMPBSA_EXE

        return engine_path

    # -----------------------------
    #           MMPBSA.in
    #
    # GB and PB calculation
    # &general
    #   startframe=1, interval=100,
    #   verbose=1, keep_files=0,
    # /
    # &gb
    #   igb=5, saltcon=0.150,
    # /
    # &pb
    #   istrng=0.15, fillratio=4.0
    # /
    CONF_IN = {
        "startframe": 1,
        "endframe": None,
        "interval": 1,
        "verbose": 1,
        "keep_files": 0,
        "if_gb": 1,
        "igb": 5,
        "saltcon": 0.15,
        "if_pb": 1,
        "istrng": 0.15,
        "fillratio": 4.0,
    }

    @classmethod
    def build_MMPBSA_in(cls, out_path=""):
        """
        build MMPBSA.in in out_path
        """
        if out_path == "":
            mkdir("./tmp")
            out_path = "./tmp/MMPBSA.in"

        # make lines
        frame_line = "  "
        for i in ("startframe", "endframe", "interval"):
            if cls.conf_in[i] != None:
                frame_line = frame_line + i + "=" + str(cls.conf_in[i]) + ", "
        output_line = (
            "  verbose="
            + str(cls.conf_in["verbose"])
            + ", keep_files="
            + str(cls.conf_in["keep_files"])
            + ","
        )
        gb_line = (
            "  igb="
            + str(cls.conf_in["igb"])
            + ", saltcon="
            + str(cls.conf_in["saltcon"])
            + ","
        )
        pb_line = (
            "  istrng="
            + str(cls.conf_in["istrng"])
            + ", fillratio="
            + str(cls.conf_in["fillratio"])
        )

        with open(out_path, "w") as of:
            print("GB and PB calculation", end=line_feed, file=of)
            print("&general", end=line_feed, file=of)
            print(frame_line, end=line_feed, file=of)
            print(output_line, end=line_feed, file=of)
            print("/", end=line_feed, file=of)
            print("&gb", end=line_feed, file=of)
            print(gb_line, end=line_feed, file=of)
            print("/", end=line_feed, file=of)
            print("&pb", end=line_feed, file=of)
            print(pb_line, end=line_feed, file=of)
            print("/", end=line_feed, file=of)

        return out_path

    def __getitem__(self, key: str) -> Any:
        return getattr(self, key)

    def __setitem__(self, key: str, value: Any) -> None:
        is_path = False
        if key in {"HOME", "EXE"}:
            is_path = True

        setattr(self, key, value)
        if is_path:
            self._parent.update_paths()
