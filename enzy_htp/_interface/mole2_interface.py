
#TODO(CJ): all of this documentation

from typing import List, Tuple

from .base_interface import BaseInterface



from enzy_htp._config.mole2_config import Mole2Config, default_mole2_config

class Mole2Cavity:
    pass

class Mole2Interface(BaseInterface):

    
    def __init__(self, parent, config: Mole2Config = None) -> None:
        """Simplistic constructor that optionally takes an AlphaFillConfig object as its only argument.
        Calls parent constructor.
        """
        super().__init__(parent, config, default_mole2_config)


    def _write_xml_input(self,  molfile:str,
                                work_dir:str,
                                non_active_parts:List[Tuple[str,int]], 
                                probe:float, 
                                inner:float, 
                                ) -> str:
        content:List[str] = [
            "<?xml version=\"1.0\" encoding=\"UTF-8\"?>", 
            "<Tunnels>",
           f"\t<WorkingDirectory>{work_dir}/</WorkingDirectory>",
           f"\t<Input>{molfile}</Input>\n"]

        if non_active_parts is not None:
            content.append("\t<NonActiveParts>")
            for (chain, rnum) in non_active_parts:
                content.append(f"\t\t<Residue Chain=\"{chain}\" SequenceNumber=\"{rnum}\" />")
        
            content.append("\t</NonActiveParts>")

        content.extend([
            "\t<Params>",
           f"\t\t<Cavity ProbeRadius=\"{probe}\" InteriorThreshold=\"{inner}\" />",
            "\t</Params>",
            "\t<Export>",
            "\t\t<Formats Mesh=\"1\" >",
            "\t\t<Mesh Density=\"0.5\" >",
            "\t<Export>",
            "</Tunnels>",
        ])

        outfile:str = f"{work_dir}/__mole2_input.xml"

        fs.write_lines(outfile, content)

        return outfile

    def identify_cavities(self, molfile:str, 
                                non_active_parts:List[Tuple[str,int]]=None, 
                                probe:float=None, 
                                inner:float=None, 
                                work_dir:str=None,
                                use_mono:bool=True
                                ) -> List[Mole2Cavity]:
        """TODO(CJ)"""

        if probe is None:
            probe = self.config_.PROBE

        if inner is None:
            inner = self.config_.INNER

        if work_dir is None:
            work_dir = config['system.SCRATCH_DIR']


        xml_file:str = self._write_xml_input(self, molfile, work_dir, non_active_parts, probe, inner )

        if use_mono:
            self.env_manager_.run_command("mono", self.config_.MOLE2, [xml_file])
        else:
            self.env_manager_.run_command(self.config_.MOLE2, [xml_file])

        
            

