from numpy import *
import waLBerla
D2Q9 = { 'C'  : 0,
         'N'  : 1,
         'S'  : 2,
         'W'  : 3,
         'E'  : 4,
         'NW' : 5,
         'NE' : 6,
         'SW' : 7,
         'SE' : 8,}


def getCollide(omega):
    def collide(block):
        wal_f = block["pdfs"] # layout: (x,y,z,nr)
        f = waLBerla.field.toArray( wal_f )
        (nx,ny,nz,dim) = f.shape
        rho = waLBerla.field.toArray( block['rho'] ).squeeze()
        u = waLBerla.field.toArray( block['vel']).squeeze()

        fpost = cumulantCollide_min(f, rho, u, omega)

        fpost = reshape(fpost,(nx,ny,nz,dim))
        copyto( asarray( block["pdfs"].buffer() ), fpost )
    return collide


def cumulantCollide_min(fin, rho, u, omega):
    (K_11, K_20, K_02) = centralMomentsFromDistributions_min ( u, fin) # the ones we need are equal

    K_11_p = (1-omega)*K_11

    K_20_p = rho/3. + 0.5*(1-omega)*(K_20 - K_02)
    K_02_p = rho/3. - 0.5*(1-omega)*(K_20 - K_02)

    return distributionsFromNormalizedCumulants_min (rho, u, K_11_p, K_20_p, K_02_p)


def centralMomentsFromDistributions_min (u, fin):
    ux = u[:, :, 0]
    uy = u[:, :, 1]
    (nx,ny) = ux.shape
    f = zeros((3,3,nx,ny))

    f[1,1] = fin[:, :, 0, D2Q9['C']]
    f[2,1] = fin[:, :, 0, D2Q9['E']]
    f[1,2] = fin[:, :, 0, D2Q9['N']]
    f[0,1] = fin[:, :, 0, D2Q9['W']]
    f[1,0] = fin[:, :, 0, D2Q9['S']]
    f[2,2] = fin[:, :, 0, D2Q9['NE']]
    f[0,2] = fin[:, :, 0, D2Q9['NW']]
    f[0,0] = fin[:, :, 0, D2Q9['SW']]
    f[2,0] = fin[:, :, 0, D2Q9['SE']]

    c_0_0 = f[0,0]              + f[0,1]           + f[0,2]
    c_0_1 = f[0,0]*(-1 - uy)    - f[0,1]*uy        + f[0,2]*(1 - uy)
    c_0_2 = f[0,0]*(-1 - uy)**2 + f[0,1]*uy**2     + f[0,2]*(1 - uy)**2
    c_1_0 = f[1,0]              + f[1,1]           + f[1,2]
    c_1_1 = f[1,0]*(-1 - uy)    - f[1,1]*uy        + f[1,2]*(1 - uy)
    c_1_2 = f[1,0]*(-1 - uy)**2 + f[1,1]*uy**2     + f[1,2]*(1 - uy)**2
    c_2_0 = f[2,0]              + f[2,1]           + f[2,2]
    c_2_1 = f[2,0]*(-1 - uy)    - f[2,1]*uy        + f[2,2]*(1 - uy)
    c_2_2 = f[2,0]*(-1 - uy)**2 + f[2,1]*uy**2     + f[2,2]*(1 - uy)**2

    c_20 = c_0_0*(-1 - ux)**2 + c_1_0*ux**2     + c_2_0*(1 - ux)**2
    c_11 = c_0_1*(-1 - ux)    - c_1_1*ux        + c_2_1*(1 - ux)
    c_02 = c_0_2              + c_1_2           + c_2_2

    return (c_11, c_20, c_02)


def distributionsFromNormalizedCumulants_min(rho, u, *normalizedCumulants):
    centralMoments = centralMomentsFromNormalizedCumulants_min (rho, *normalizedCumulants)

    return distributionsFromCentralMoments_min (u, *centralMoments)


def centralMomentsFromNormalizedCumulants_min (rho, K_11, K_20, K_02):
    # Transformation to central moments
    c_22 = 2*K_11*K_11/rho + K_20*K_02/rho

    return (rho, K_11, K_20, K_02, c_22)


def distributionsFromCentralMoments_min (u, c_00, c_11, c_20, c_02, c_22):
    ux = u[:, :, 0]
    uy = u[:, :, 1]
    (nx,ny) = ux.shape
    f = zeros((nx,ny,1,9))
    c_0_0 = (c_00*(ux**2 - ux)  + c_20) * 0.5
    c_1_0 =  c_00*(1 - ux**2)   - c_20
    c_2_0 = (c_00*(ux**2 + ux)  + c_20) * 0.5
    c_0_1 = ( c_11*(2*ux - 1) ) * 0.5
    c_1_1 =  - c_11*2*ux
    c_2_1 = ( c_11*(2*ux + 1) ) * 0.5
    c_0_2 = (c_02*(ux**2 - ux) + c_22) * 0.5
    c_1_2 =  c_02*(1 - ux**2)  - c_22
    c_2_2 = (c_02*(ux**2 + ux) + c_22) * 0.5
    f[:, :, 0, D2Q9['SW']] = (c_0_0*(uy**2 - uy) + c_0_1*(2*uy - 1) + c_0_2) * 0.5 # (0, 0)
    f[:, :, 0, D2Q9['W']]  =  c_0_0*(1 - uy**2)  - c_0_1*2*uy       - c_0_2        # (0, 1)
    f[:, :, 0, D2Q9['NW']] = (c_0_0*(uy**2 + uy) + c_0_1*(2*uy + 1) + c_0_2) * 0.5 # (0, 2)
    f[:, :, 0, D2Q9['S']]  = (c_1_0*(uy**2 - uy) + c_1_1*(2*uy - 1) + c_1_2) * 0.5 # (1, 0)
    f[:, :, 0, D2Q9['C']]  =  c_1_0*(1 - uy**2)  - c_1_1*2*uy       - c_1_2        # (1, 1)
    f[:, :, 0, D2Q9['N']]  = (c_1_0*(uy**2 + uy) + c_1_1*(2*uy + 1) + c_1_2) * 0.5 # (1, 2)
    f[:, :, 0, D2Q9['SE']] = (c_2_0*(uy**2 - uy) + c_2_1*(2*uy - 1) + c_2_2) * 0.5 # (2, 0)
    f[:, :, 0, D2Q9['E']]  =  c_2_0*(1 - uy**2)  - c_2_1*2*uy       - c_2_2        # (2, 1)
    f[:, :, 0, D2Q9['NE']] = (c_2_0*(uy**2 + uy) + c_2_1*(2*uy + 1) + c_2_2) * 0.5 # (2, 2)

    return f

def getDensity(f,nx,ny):
    rho = f[:,:,0,0] + f[:,:,0,1] + f[:,:,0,2] + f[:,:,0,3] + f[:,:,0,4] + f[:,:,0,5] + f[:,:,0,6] + f[:,:,0,7] + f[:,:,0,8]
    return rho
