#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include "contact_map.h"

// REFERENCE: J. Tsai, R. Taylor, C. Chothia, and M. Gerstein, J. Mol. Biol 290:290 (1999)
// REFERENCE: https://aip.scitation.org/doi/suppl/10.1063/1.4929599/suppl_file/sm.pdf

// ALL ALPHABET RELATED TO PROTEINS
typedef enum {
	ACE=1,ALA,ARG,ASN,ASP,CLA,CYS,GLN,GLU,GLY,HIS,HYP,ILE,L4Y,L5Y,LEU,LY2,LY3,LYS,LYX,MET,NME,PHE,PRO,SER,THR,TRP,TYR,VAL,AGS,LGX,LZS,LZD,
	CYM,CYX,HID,HIE,HIP,
	C,CA,CB,CD,CD1,CD2,CE,CE1,CE2,CE3,CG,CG1,CG2,CH2,CZ,CZ2,CZ3,CH3,C10,C11,C12,C13,C14,C15,C16,C17,C18,C19,C20,C28,C30,C31,
	N,ND1,ND2,NE,NE1,NE2,NH1,NH2,NZ,NZ1,NZ2,
	O,OD1,OD2,OE1,OE2,OG,OG1,OH,OE,OD,O11,O18,O19,OXT,
	SD,SG,
} protein_enum;

// PROTEIN LOOKUP TABLE
typedef struct {
	char *key;
	protein_enum val;
} protein_str;

static protein_str plt[] = {
	{"ACE",ACE},{"ALA",ALA},{"ARG",ARG},{"ASN",ASN},{"ASP",ASP},{"CLA",CLA},{"CYS",CYS},
	{"GLN",GLN},{"GLU",GLU},{"GLY",GLY},{"HIS",HIS},{"HYP",HYP},{"ILE",ILE},{"L4Y",L4Y},
        {"L5Y",L5Y},{"LEU",LEU},{"LY2",LY2},{"LY3",LY3},{"LYS",LYS},{"LYX",LYX},{"MET",MET},
        {"NME",NME},{"PHE",PHE},{"PRO",PRO},{"SER",SER},{"THR",THR},{"TRP",TRP},{"TYR",TYR},
        {"VAL",VAL},{"AGS",AGS},{"LGX",LGX},{"LZS",LZS},{"LZD",LZD},{"CYM",CYM},{"CYX",CYX},{"HID",HID},{"HIE",HIE},{"HIP",HIP},
	{"C"  ,  C},{"CA" , CA},{"CB" , CB},{"CD" , CD},{"CD1",CD1},{"CD2",CD2},{"CE" , CE},
	{"CE1",CE1},{"CE2",CE2},{"CE3",CE3},{"CG" , CG},{"CG1",CG1},{"CG2",CG2},{"CH2",CH2},
	{"CZ" , CZ},{"CZ2",CZ2},{"CZ3",CZ3},{"CH3",CH3},{"C10",C10},{"C11",C11},{"C12",C12},{"C13",C13},{"C14",C14},
	{"C15",C15},{"C16",C16},{"C17",C17},{"C18",C18},{"C19",C19},{"C20",C20},
	{"C28",C28},{"C30",C30},{"C31",C31},
	{"N"  ,  N},{"ND1",ND1},{"ND2",ND2},{"NE" , NE},{"NE1",NE1},{"NE2",NE2},{"NH1",NH1},
	{"NH2",NH2},{"NZ" , NZ},{"NZ1",NZ1},{"NZ2",NZ2},
	{"O"  ,  O},{"OD1",OD1},{"OD2",OD2},{"OE1",OE1},{"OE2",OE2},{"OG" , OG},{"OG1",OG1},
	{"OH" , OH},{"OE",  OE},{"OD",  OD},{"O11",O11},{"O18",O18},{"O19",O19},{"OXT",OXT},
	{"SD" , SD},{"SG" , SG}
};

static int keyfromstring(char key[]) {
	int r = 0;
	for (int k=0;k<sizeof(plt)/sizeof(protein_str);k++)
		if ( (strlen(plt[k].key) == strlen(key)) && (strncmp(plt[k].key,key,strlen(key)) == 0) ) {
			r=plt[k].val;
			break;
		}
	//if (r<=0) {
	//	printf("Key '%s' not found on protein map!\n",key);
	//};
	return r;
}

// THIS SUBROUTINE ASSIGNS VAN DER WAALS RADII
extern bool protein_map(atom_pdb_str *atom, atomaux_str *vdw) {

	vdw->keyatom = keyfromstring(atom->name);
	vdw->keyresidue = keyfromstring(atom->resName);

	float vrad = 0.0;
	int atype = 0;

	switch(keyfromstring(atom->resName)) {
		case ACE:
			vdw->nb=3;
			switch(keyfromstring(atom->name)) {
				case C  : vrad=1.61; atype=6; break;
				case O  : vrad=1.42; atype=2; break;
				case CH3: vrad=1.88; atype=4; break;
				default : vrad=0.00; atype=0; printf("UNMAPPED ATOM %s %s\n",atom->resName,atom->name);
			} break;
		case ALA:
			vdw->nb=6;
			switch(keyfromstring(atom->name)) {
				case N  : vrad=1.64; atype=3; break;
				case CA : vrad=1.88; atype=7; break;
				case C  : vrad=1.61; atype=6; break;
				case O  : vrad=1.42; atype=2; break;
				case CB : vrad=1.88; atype=4; break;
				case OXT: vrad=1.42; atype=2; break;
				default : vrad=0.00; atype=0; printf("UNMAPPED ATOM %s %s\n",atom->resName,atom->name);
			} break;
		case ARG:
			vdw->nb=11;
			switch(keyfromstring(atom->name)) {
				case N  : vrad=1.64; atype=3; break;
				case CA : vrad=1.88; atype=7; break;
				case C  : vrad=1.61; atype=6; break;
				case O  : vrad=1.42; atype=2; break;
				case CB : vrad=1.88; atype=4; break;
				case CG : vrad=1.88; atype=4; break;
				case CD : vrad=1.88; atype=7; break;
				case NE : vrad=1.64; atype=3; break;
				case CZ : vrad=1.61; atype=6; break;
				case NH1: vrad=1.64; atype=3; break;
				case NH2: vrad=1.64; atype=3; break;
				default : vrad=0.00; atype=0; printf("UNMAPPED ATOM %s %s\n",atom->resName,atom->name);
			} break;
		case AGS:
			vdw->nb=20;
			switch(keyfromstring(atom->name)) {
				case N  : vrad=1.64; atype=3; break;
				case CA : vrad=1.88; atype=7; break;
				case C  : vrad=1.61; atype=6; break;
				case O  : vrad=1.42; atype=2; break;
				case CB : vrad=1.88; atype=4; break;
				case CG : vrad=1.88; atype=4; break;
				case CD : vrad=1.88; atype=7; break;
				case NE : vrad=1.64; atype=3; break;
				case CZ : vrad=1.61; atype=6; break;
				case NH1: vrad=1.64; atype=2; break;
				case NH2: vrad=1.64; atype=2; break;
				case NZ : vrad=1.64; atype=2; break;
				/* Glucosepane is non-aromatic: classify its saturated ring
				 * carbons by their directly bonded heteroatoms, not as type 5. */
				case C15: vrad=1.88; atype=6; break;
				case C16: vrad=1.61; atype=6; break;
				case C17: vrad=1.88; atype=7; break;
				case C18: vrad=1.88; atype=6; break;
				case O18: vrad=1.46; atype=1; break;
				case C19: vrad=1.88; atype=6; break;
				case O19: vrad=1.46; atype=1; break;
				case C20: vrad=1.88; atype=4; break;
				default : vrad=0.00; atype=0; printf("UNMAPPED ATOM %s %s\n",atom->resName,atom->name);
			} break;
		case ASN:
			vdw->nb=8;
			switch(keyfromstring(atom->name)) {
				case N  : vrad=1.64; atype=3; break;
				case CA : vrad=1.88; atype=7; break;
				case C  : vrad=1.61; atype=6; break;
				case O  : vrad=1.42; atype=2; break;
				case CB : vrad=1.88; atype=4; break;
				case CG : vrad=1.61; atype=6; break;
				case OD1: vrad=1.42; atype=2; break;
				case ND2: vrad=1.64; atype=3; break;
				default : vrad=0.00; atype=0; printf("UNMAPPED ATOM %s %s\n",atom->resName,atom->name);
			} break;
		case ASP:
			vdw->nb=8;
			switch(keyfromstring(atom->name)) {
				case N  : vrad=1.64; atype=3; break;
				case CA : vrad=1.88; atype=7; break;
				case C  : vrad=1.61; atype=6; break;
				case O  : vrad=1.42; atype=2; break;
				case CB : vrad=1.88; atype=4; break;
				case CG : vrad=1.61; atype=6; break;
				case OD1: vrad=1.46; atype=2; break;
				case OD2: vrad=1.42; atype=2; break;
				default : vrad=0.00; atype=0; printf("UNMAPPED ATOM %s %s\n",atom->resName,atom->name);
			} break;
		case CLA:
                        vdw->nb=6;
                        switch(keyfromstring(atom->name)) {
                                case N  : vrad=1.64; atype=3; break;
                                case CA : vrad=1.88; atype=7; break;
                                case C  : vrad=1.61; atype=6; break;
                                case O  : vrad=1.42; atype=2; break;
                                case CB : vrad=1.88; atype=4; break;
                                case OXT: vrad=1.42; atype=2; break;
                                default : vrad=0.00; atype=0; printf("UNMAPPED ATOM %s %s\n",atom->resName,atom->name);
                        } break;
		case CYM:
		case CYX:
		case CYS:
			vdw->nb=6;
			switch(keyfromstring(atom->name)) {
				case N  : vrad=1.64; atype=3; break;
				case CA : vrad=1.88; atype=7; break;
				case C  : vrad=1.61; atype=6; break;
				case O  : vrad=1.42; atype=2; break;
				case CB : vrad=1.88; atype=4; break;
				case SG : vrad=1.77; atype=6; break;
				default : vrad=0.00; atype=0; printf("UNMAPPED ATOM %s %s\n",atom->resName,atom->name);
			} break;
		case GLN:
			vdw->nb=9;
			switch(keyfromstring(atom->name)) {
				case N  : vrad=1.64; atype=3; break;
				case CA : vrad=1.88; atype=7; break;
				case C  : vrad=1.61; atype=6; break;
				case O  : vrad=1.42; atype=2; break;
				case CB : vrad=1.88; atype=4; break;
				case CG : vrad=1.88; atype=4; break;
				case CD : vrad=1.61; atype=6; break;
				case OE1: vrad=1.42; atype=2; break;
				case NE2: vrad=1.64; atype=3; break;
				default : vrad=0.00; atype=0; printf("UNMAPPED ATOM %s %s\n",atom->resName,atom->name);
			} break;
		case GLU:
			vdw->nb=9;
			switch(keyfromstring(atom->name)) {
				case N  : vrad=1.64; atype=3; break;
				case CA : vrad=1.88; atype=7; break;
				case C  : vrad=1.61; atype=6; break;
				case O  : vrad=1.42; atype=2; break;
				case CB : vrad=1.88; atype=4; break;
				case CG : vrad=1.88; atype=4; break;
				case CD : vrad=1.61; atype=6; break;
				case OE1: vrad=1.46; atype=2; break;
				case OE2: vrad=1.42; atype=2; break;
				default : vrad=0.00; atype=0; printf("UNMAPPED ATOM %s %s\n",atom->resName,atom->name);
			} break;
		case GLY:
			vdw->nb=4;
			switch(keyfromstring(atom->name)) {
				case N  : vrad=1.64; atype=3; break;
				case CA : vrad=1.88; atype=6; break;
				case C  : vrad=1.61; atype=6; break;
				case O  : vrad=1.42; atype=2; break;
				default : vrad=0.00; atype=0; printf("UNMAPPED ATOM %s %s\n",atom->resName,atom->name);
			} break;
		case HIE:
		case HIP:
		case HYP:
			vdw->nb=8;
			switch(keyfromstring(atom->name)) {
				case N  : vrad=1.64; atype=6; break;
				case CA : vrad=1.88; atype=4; break;
				case C  : vrad=1.61; atype=6; break;
				case O  : vrad=1.42; atype=2; break;
				case CB : vrad=1.88; atype=4; break;
				case CG : vrad=1.88; atype=4; break;
				case CD2: vrad=1.76; atype=6; break;
				case OD1: vrad=1.46; atype=1; break;
				default : vrad=0.00; atype=0; printf("UNMAPPED ATOM %s %s\n",atom->resName,atom->name);
			} break;
		case HIS:
			vdw->nb=10;
			switch(keyfromstring(atom->name)) {
				case N  : vrad=1.64; atype=3; break;
				case CA : vrad=1.88; atype=7; break;
				case C  : vrad=1.61; atype=6; break;
				case O  : vrad=1.42; atype=2; break;
				case CB : vrad=1.88; atype=4; break;
				case CG : vrad=1.61; atype=5; break;
				case ND1: vrad=1.64; atype=1; break;
				case CD2: vrad=1.76; atype=5; break;
				case CE1: vrad=1.76; atype=5; break;
				case NE2: vrad=1.64; atype=1; break;
				default : vrad=0.00; atype=0; printf("UNMAPPED ATOM %s %s\n",atom->resName,atom->name);
			} break;
		case ILE:
			vdw->nb=8;
			switch(keyfromstring(atom->name)) {
				case N  : vrad=1.64; atype=3; break;
				case CA : vrad=1.88; atype=7; break;
				case C  : vrad=1.61; atype=6; break;
				case O  : vrad=1.42; atype=2; break;
				case CB : vrad=1.88; atype=4; break;
				case CG1: vrad=1.88; atype=4; break;
				case CG2: vrad=1.88; atype=4; break;
				case CD1: vrad=1.88; atype=4; break;
				case  CD: vrad=1.88; atype=4; break;
				default : vrad=0.00; atype=0; printf("UNMAPPED ATOM %s %s\n",atom->resName,atom->name);
			} break;
		case L4Y:
			vdw->nb=9;
			switch(keyfromstring(atom->name)) {
				case N  : vrad=1.64; atype=3; break;
				case CA : vrad=1.88; atype=7; break;
				case C  : vrad=1.61; atype=6; break;
				case O  : vrad=1.42; atype=2; break;
				case CB : vrad=1.88; atype=4; break;
				case CG : vrad=1.88; atype=4; break;
				case CD : vrad=1.88; atype=4; break;
				case OE : vrad=1.42; atype=2; break;
				case CE : vrad=1.61; atype=6; break;
				default : vrad=0.00; atype=0; printf("UNMAPPED ATOM %s %s\n",atom->resName,atom->name);
			} break;
		case L5Y:
			vdw->nb=10;
			switch(keyfromstring(atom->name)) {
				case N  : vrad=1.64; atype=3; break;
				case CA : vrad=1.88; atype=7; break;
				case C  : vrad=1.61; atype=6; break;
				case O  : vrad=1.42; atype=2; break;
				case CB : vrad=1.88; atype=4; break;
				case CG : vrad=1.88; atype=4; break;
				case CD : vrad=1.76; atype=6; break;
				case OD : vrad=1.46; atype=1; break;
				case CE : vrad=1.76; atype=6; break;
				case NZ : vrad=1.64; atype=1; break;
				default : vrad=0.00; atype=0; printf("UNMAPPED ATOM %s %s\n",atom->resName,atom->name);
			} break;
		case LEU:
			vdw->nb=8;
			switch(keyfromstring(atom->name)) {
				case N  : vrad=1.64; atype=3; break;
				case CA : vrad=1.88; atype=7; break;
				case C  : vrad=1.61; atype=6; break;
				case O  : vrad=1.42; atype=2; break;
				case CB : vrad=1.88; atype=4; break;
				case CG : vrad=1.88; atype=4; break;
				case CD1: vrad=1.88; atype=4; break;
				case CD2: vrad=1.88; atype=4; break;
				default : vrad=0.00; atype=0; printf("UNMAPPED ATOM %s %s\n",atom->resName,atom->name);
			} break;
		case LY2:
			vdw->nb=5;
			switch(keyfromstring(atom->name)) {
				case N  : vrad=1.64; atype=3; break;
				case CA : vrad=1.88; atype=7; break;
				case C  : vrad=1.61; atype=6; break;
				case O  : vrad=1.42; atype=2; break;
				case CB : vrad=1.88; atype=4; break;
				default : vrad=0.00; atype=0; printf("UNMAPPED ATOM %s %s\n",atom->resName,atom->name);
			} break;
		case LY3:
			vdw->nb=6;
			switch(keyfromstring(atom->name)) {
				case N  : vrad=1.64; atype=3; break;
				case CA : vrad=1.88; atype=7; break;
				case C  : vrad=1.61; atype=6; break;
				case O  : vrad=1.42; atype=2; break;
				case CB : vrad=1.88; atype=4; break;
				case CG : vrad=1.88; atype=4; break;
				default : vrad=0.00; atype=0; printf("UNMAPPED ATOM %s %s\n",atom->resName,atom->name);
			} break;
		case LYX:
			vdw->nb=16;
			switch(keyfromstring(atom->name)) {
				case N  : vrad=1.64; atype=3; break;
				case CA : vrad=1.88; atype=7; break;
				case C  : vrad=1.61; atype=6; break;
				case O  : vrad=1.42; atype=2; break;
				case CB : vrad=1.88; atype=4; break;
				case CG : vrad=1.88; atype=4; break;
				case CD : vrad=1.76; atype=6; break;
				case OD : vrad=1.46; atype=1; break;
				case CE : vrad=1.76; atype=7; break;
				case NZ : vrad=1.46; atype=1; break;
				case C10: vrad=1.76; atype=5; break;
				case C11: vrad=1.76; atype=5; break;
				case O11: vrad=1.46; atype=1; break;
				case C12: vrad=1.76; atype=5; break;
				case C13: vrad=1.76; atype=5; break;
				case C14: vrad=1.76; atype=5; break;
				default : vrad=0.00; atype=0; printf("UNMAPPED ATOM %s %s\n",atom->resName,atom->name);
			} break;
		case LYS:
			vdw->nb=9;
			switch(keyfromstring(atom->name)) {
				case N  : vrad=1.64; atype=3; break;
				case CA : vrad=1.88; atype=7; break;
				case C  : vrad=1.61; atype=6; break;
				case O  : vrad=1.42; atype=2; break;
				case CB : vrad=1.88; atype=4; break;
				case CG : vrad=1.88; atype=4; break;
				case CD : vrad=1.88; atype=4; break;
				case CE : vrad=1.88; atype=7; break;
				case NZ : vrad=1.64; atype=3; break;
				default : vrad=0.00; atype=0; printf("UNMAPPED ATOM %s %s\n",atom->resName,atom->name);
			} break;
		case LGX:
			vdw->nb=8;
			switch(keyfromstring(atom->name)) {
				case N  : vrad=1.64; atype=3; break;
				case CA : vrad=1.88; atype=7; break;
				case C  : vrad=1.61; atype=6; break;
				case O  : vrad=1.42; atype=2; break;
				case CB : vrad=1.88; atype=4; break;
				case CG : vrad=1.88; atype=4; break;
				case CD : vrad=1.88; atype=4; break;
				case CE : vrad=1.88; atype=7; break;
				default : vrad=0.00; atype=0; printf("UNMAPPED ATOM %s %s\n",atom->resName,atom->name);
			} break;
		case LZS:
			vdw->nb=14;
			switch(keyfromstring(atom->name)) {
				case N  : vrad=1.64; atype=3; break;
				case CA : vrad=1.88; atype=7; break;
				case C  : vrad=1.61; atype=6; break;
				case O  : vrad=1.42; atype=2; break;
				case CB : vrad=1.88; atype=4; break;
				case CG : vrad=1.88; atype=4; break;
				case CD : vrad=1.88; atype=4; break;
				case CE : vrad=1.88; atype=7; break;
				case NZ1: vrad=1.64; atype=3; break;
				case NZ2: vrad=1.64; atype=3; break;
				/* MOLD contains an aromatic imidazolium ring. */
				case C28: vrad=1.76; atype=5; break;
				case C30: vrad=1.76; atype=5; break;
				case C31: vrad=1.76; atype=5; break;
				case CH3: vrad=1.88; atype=4; break;
				default : vrad=0.00; atype=0; printf("UNMAPPED ATOM %s %s\n",atom->resName,atom->name);
			} break;
		case LZD:
			vdw->nb=8;
			switch(keyfromstring(atom->name)) {
				case N  : vrad=1.64; atype=3; break;
				case CA : vrad=1.88; atype=7; break;
				case C  : vrad=1.61; atype=6; break;
				case O  : vrad=1.42; atype=2; break;
				case CB : vrad=1.88; atype=4; break;
				case CG : vrad=1.88; atype=4; break;
				case CD : vrad=1.88; atype=4; break;
				case CE : vrad=1.88; atype=7; break;
				default : vrad=0.00; atype=0; printf("UNMAPPED ATOM %s %s\n",atom->resName,atom->name);
			} break;
		case MET:
			vdw->nb=8;
			switch(keyfromstring(atom->name)) {
				case N  : vrad=1.64; atype=3; break;
				case CA : vrad=1.88; atype=7; break;
				case C  : vrad=1.61; atype=6; break;
				case O  : vrad=1.42; atype=2; break;
				case CB : vrad=1.88; atype=4; break;
				case CG : vrad=1.88; atype=4; break;
				case SD : vrad=1.77; atype=8; break;
				case CE : vrad=1.88; atype=4; break;
				default : vrad=0.00; atype=0; printf("UNMAPPED ATOM %s %s\n",atom->resName,atom->name);
			} break;
		case NME:
			vdw->nb=2;
			switch(keyfromstring(atom->name)) {
				case N  : vrad=1.64; atype=3; break;
				case CH3: vrad=1.88; atype=4; break;
				default : vrad=0.00; atype=0; printf("UNMAPPED ATOM %s %s\n",atom->resName,atom->name);
			} break;
		case PHE:
			vdw->nb=11;
			switch(keyfromstring(atom->name)) {
				case N  : vrad=1.64; atype=3; break;
				case CA : vrad=1.88; atype=7; break;
				case C  : vrad=1.61; atype=6; break;
				case O  : vrad=1.42; atype=2; break;
				case CB : vrad=1.88; atype=4; break;
				case CG : vrad=1.88; atype=5; break;
				case CD1: vrad=1.61; atype=5; break;
				case CD2: vrad=1.76; atype=5; break;
				case CE1: vrad=1.76; atype=5; break;
				case CE2: vrad=1.76; atype=5; break;
				case CZ : vrad=1.76; atype=5; break;
				default : vrad=0.00; atype=0; printf("UNMAPPED ATOM %s %s\n",atom->resName,atom->name);
			} break;
		case PRO:
			vdw->nb=7;
			switch(keyfromstring(atom->name)) {
				case N  : vrad=1.64; atype=6; break;
				case CA : vrad=1.88; atype=4; break;
				case C  : vrad=1.61; atype=6; break;
				case O  : vrad=1.42; atype=2; break;
				case CB : vrad=1.88; atype=4; break;
				case CG : vrad=1.88; atype=4; break;
				case CD : vrad=1.88; atype=4; break;
				default : vrad=0.00; atype=0; printf("UNMAPPED ATOM %s %s\n",atom->resName,atom->name);
			} break;
		case SER:
			vdw->nb=6;
			switch(keyfromstring(atom->name)) {
				case N  : vrad=1.64; atype=3; break;
				case CA : vrad=1.88; atype=7; break;
				case C  : vrad=1.61; atype=6; break;
				case O  : vrad=1.42; atype=2; break;
				case CB : vrad=1.88; atype=6; break;
				case OG : vrad=1.46; atype=1; break;
				default : vrad=0.00; atype=0; printf("UNMAPPED ATOM %s %s\n",atom->resName,atom->name);
			} break;
		case THR:
			vdw->nb=7;
			switch(keyfromstring(atom->name)) {
				case N  : vrad=1.64; atype=3; break;
				case CA : vrad=1.88; atype=7; break;
				case C  : vrad=1.61; atype=6; break;
				case O  : vrad=1.42; atype=2; break;
				case CB : vrad=1.88; atype=6; break;
				case OG1: vrad=1.46; atype=1; break;
				case CG2: vrad=1.88; atype=4; break;
				default : vrad=0.00; atype=0; printf("UNMAPPED ATOM %s %s\n",atom->resName,atom->name);
			} break;
		case TRP:
			vdw->nb=14;
			switch(keyfromstring(atom->name)) {
				case N  : vrad=1.64; atype=3; break;
				case CA : vrad=1.88; atype=7; break;
				case C  : vrad=1.61; atype=6; break;
				case O  : vrad=1.42; atype=2; break;
				case CB : vrad=1.88; atype=4; break;
				case CG : vrad=1.61; atype=5; break;
				case CD1: vrad=1.76; atype=5; break;
				case CD2: vrad=1.61; atype=5; break;
				case NE1: vrad=1.64; atype=3; break;
				case CE2: vrad=1.61; atype=5; break;
				case CE3: vrad=1.76; atype=5; break;
				case CZ2: vrad=1.76; atype=5; break;
				case CZ3: vrad=1.76; atype=5; break;
				case CH2: vrad=1.76; atype=5; break;
				default : vrad=0.00; atype=0; printf("UNMAPPED ATOM %s %s\n",atom->resName,atom->name);
			} break;
		case TYR:
			vdw->nb=12;
			switch(keyfromstring(atom->name)) {
				case N  : vrad=1.64; atype=3; break;
				case CA : vrad=1.88; atype=7; break;
				case C  : vrad=1.61; atype=6; break;
				case O  : vrad=1.42; atype=2; break;
				case CB : vrad=1.88; atype=4; break;
				case CG : vrad=1.61; atype=5; break;
				case CD1: vrad=1.76; atype=5; break;
				case CD2: vrad=1.76; atype=5; break;
				case CE1: vrad=1.76; atype=5; break;
				case CE2: vrad=1.76; atype=5; break;
				case CZ : vrad=1.61; atype=5; break;
				case OH : vrad=1.46; atype=1; break;
				default : vrad=0.00; atype=0; printf("UNMAPPED ATOM %s %s\n",atom->resName,atom->name);
			} break;
		case VAL:
			vdw->nb=7;
			switch(keyfromstring(atom->name)) {
				case N  : vrad=1.64; atype=3; break;
				case CA : vrad=1.88; atype=7; break;
				case C  : vrad=1.61; atype=6; break;
				case O  : vrad=1.42; atype=2; break;
				case CB : vrad=1.88; atype=4; break;
				case CG1: vrad=1.88; atype=4; break;
				case CG2: vrad=1.88; atype=4; break;
				default : vrad=0.00; atype=0; printf("UNMAPPED ATOM %s %s\n",atom->resName,atom->name);
			} break;
		default:
			printf("UNMAPPED RESIDUE %s %s\n",atom->resName,atom->name);
	} // switch RESIDUE
	vdw->vrad = vrad;
	vdw->atype = atype;
	return vdw->vrad == 0.0 && vdw->atype == 0 ? false : true;
}
