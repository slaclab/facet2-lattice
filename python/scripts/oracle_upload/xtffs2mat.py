#!/bin/env python3

import numpy as np

def xtffs2mat(fname):
    """
    Reads an XTFF SURVEY file and extracts various parameters.

    Parameters:
    fname (str): The filename of the XTFF SURVEY file.

    Returns:
    tuple: Contains the following elements:
        tt (str): Run title
        K (list): Element keyword
        N (list): Element name
        L (np.array): Element length
        P (np.array): Element parameter
        A (np.array): Aperture
        T (list): Engineering type
        E (np.array): Energy
        FDN (list): NLC Formal Device Name
        coor (np.array): Survey coordinates (X,Y,Z,yaw,pitch,roll)
        S (np.array): suml
    """

    # Open the XTFF file
    try:
        with open(fname, 'r') as fid:
            # Read in the header ... check that XTFF file is a SURVEY file
            line = fid.readline()
            xtff = line[8:16]
            if xtff.strip() != "SURVEY":
                raise ValueError(f"Unexpected XTFF type ({xtff}) encountered ... abort")

            # Read in the run title
            tt = fid.readline().strip()

            # Initialize lists and arrays
            K, N, L, P, A, T, E, FDN, coor, S = [], [], [], [], [], [], [], [], [], []

            # Read in the INITIAL data
            for _ in range(4):
                line = fid.readline()
                if _ == 0:
                    K.append(line[0:4].strip())
                    N.append(line[4:20].strip())
                    L.append(float(line[20:32].strip()))
                    p = [float(line[32:48].strip()), float(line[48:64].strip()), float(line[64:80].strip())]
                    A.append(float(line[80:96].strip()))
                    T.append(line[97:113].strip())
                    E.append(float(line[114:130].strip()))
                elif _ == 1:
                    line = line.ljust(105)
                    p.extend([float(line[0:16].strip()), float(line[16:32].strip()), float(line[32:48].strip()),
                              float(line[48:64].strip()), float(line[64:80].strip())])
                    FDN.append(line[81:105].strip())
                elif _ == 2:
                    x, y, z = map(float, [line[0:16].strip(), line[16:32].strip(), line[32:48].strip()])
                    S.append(float(line[48:64].strip()))
                else:
                    yaw, pitch, roll = map(float, [line[0:16].strip(), line[16:32].strip(), line[32:48].strip()])

            P.append(p)
            coor.append([x, y, z, yaw, pitch, roll])

            # Read in the data ... break at end of the file
            while True:
                line = fid.readline()
                if not line or not line[0].isalpha():
                    break

                K.append(line[0:4].strip())
                N.append(line[4:20].strip())
                L.append(float(line[20:32].strip()))
                p = [float(line[32:48].strip()), float(line[48:64].strip()), float(line[64:80].strip())]
                A.append(float(line[80:96].strip()))
                T.append(line[97:113].strip())
                E.append(float(line[114:130].strip()))

                line = fid.readline().ljust(105)
                p.extend([float(line[0:16].strip()), float(line[16:32].strip()), float(line[32:48].strip()),
                          float(line[48:64].strip()), float(line[64:80].strip())])
                FDN.append(line[81:105].strip())

                line = fid.readline()
                x, y, z = map(float, [line[0:16].strip(), line[16:32].strip(), line[32:48].strip()])
                S.append(float(line[48:64].strip()))

                line = fid.readline()
                yaw, pitch, roll = map(float, [line[0:16].strip(), line[16:32].strip(), line[32:48].strip()])

                P.append(p)
                coor.append([x, y, z, yaw, pitch, roll])

    except IOError:
        raise IOError(f"Failed to open {fname}")

    # Convert lists to numpy arrays where appropriate
    L = np.array(L)
    P = np.array(P)
    A = np.array(A)
    E = np.array(E)
    coor = np.array(coor)
    S = np.array(S)

    return tt, K, N, L, P, A, T, E, FDN, coor, S


