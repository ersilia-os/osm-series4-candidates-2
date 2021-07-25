"""Wrapper for ``mol2svg`` tool."""
import os
import subprocess
import shutil

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import cm

from rdkit import Chem

from tqdm import tqdm

class MplColorHelper:

    def __init__(self, cmap_name, start_val, stop_val):
        self.cmap_name = cmap_name
        self.cmap = plt.get_cmap(cmap_name)
        self.norm = mpl.colors.Normalize(vmin=start_val, vmax=stop_val)
        self.scalarMap = cm.ScalarMappable(norm=self.norm, cmap=self.cmap)

    def get_rgb(self, val):
        return self.scalarMap.to_rgba(val)


class Mol2svg():
    """Mol2svg class."""

    def __init__(self, output_path, cmap_name="Spectral"):
        self.exec_path = os.path.join(os.path.dirname(
            os.path.abspath(__file__)), "files")
        self.output_path = output_path
        if os.path.exists(output_path):
            shutil.rmtree(output_path)
        os.makedirs(output_path, exist_ok=True)
        self.cmap_name = cmap_name
        self.cmap = MplColorHelper(cmap_name, 0, 1)

    def mol_to_mol_file(self, mol, name):
        molObj = mol
        molObj.SetProp("_Name", name)
        imported = Chem.MolToMolBlock(molObj)
        with open("{0}/{1}.mol".format(self.output_path, name), "w") as newfile:
            newfile.write(imported)

    def mol2svg(self, mol, name, value=None):
        """Convert molecule to SVG

        Args:
            output_path:Output data path
        """
        # check input
        if value is None:
            rgb = (255, 255, 255)
        else:
            rgb = self.cmap.get_rgb(value)
            rgb = (int(rgb[0]*255), int(rgb[1]*255), int(rgb[2]*255))

        self.mol_to_mol_file(mol, name)

        rgb = "%d,%d,%d" % rgb

        cmd = '''
        {0}/mol2svg \
           --bgcolor={1} \
           --color={0}/black.conf \
           {2}/{3}.mol > \
           {2}/{3}.svg
        '''.format(
            self.exec_path,
            rgb,
            self.output_path,
            name
        )

        # run process
        process = subprocess.Popen(cmd, shell=True,
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)

    def clean_mol_files(self):
        for f in os.listdir(self.output_path):
            if f[-3:] == "mol":
                os.remove(os.path.join(self.output_path, f))

    def modify_svg(self, svg_file):
        with open(svg_file, "r") as f:
            svg = f.readlines()

        default_width  = 350
        default_height = 350

        i = 1
        l = svg[i]
        width = int(l.split('width="')[1].split('"')[0])
        height = int(l.split('height="')[1].split('"')[0])
        l = l.replace('width="%d"' % width, 'width="%d"' % default_width)
        l = l.replace('height="%d"' % height, 'height="%d"' % default_height)
        svg[i] = l

        i = 6
        l = svg[i]
        l = l.replace('y="5"', 'y="0"')
        l = l.replace('width="%d"' % width, 'width="%d"' % default_width)
        l = l.replace('height="%d"' % height, 'height="%d"' % default_height)
        svg[i] = l

        i = 8
        l = svg[i]
        if "translate" not in l:
            l = '<g transform="translate(0,0)">\n'
        x,y = l.split("translate(")[1].split(")")[0].split(",")
        x,y = int(x), int(y)
        xt = x + (default_width - width)/2
        yt = y + (default_height - height)/2
        l = l.replace('translate(%d,%d)' % (x,y), 'translate(%d,%d)' % (xt,yt))
        svg[i] = l

        with open(svg_file, "w") as f:
            f.write("".join(svg))

    def modify_svgs(self):
        for f in os.listdir(self.output_path):
            if f[-3:] == "svg":
                try:
                    self.modify_svg(os.path.join(self.output_path, f))
                except:
                    os.remove(os.path.join(self.output_path, f))

    def convert_to_png(self):
        for f in os.listdir(self.output_path):
            if f[-3:] == "svg":
                fn = f.split(".")[0]
                cmd = "inkscape -z -w 480 -h 480 {0}/{1}.svg -e {0}/{1}.png".format(self.output_path, fn)
                process = subprocess.Popen(cmd, shell=True,
                                           stdout=subprocess.PIPE,
                                           stderr=subprocess.PIPE)

    def clean_svg_files(self):
        for f in os.listdir(self.output_path):
            if f[-3:] == "svg":
                os.remove(os.path.join(self.output_path, f))

    def mols2png(self, mols, values=None):
        i = 0
        for i, mol in tqdm(enumerate(mols)):
            name = "mol%05d" % i
            if values is None:
                v = None
            else:
                v = values[i]
            self.mol2svg(mol, name, v)
        self.clean_mol_files()
        self.modify_svgs()
        self.convert_to_png()
        self.clean_svg_files()
