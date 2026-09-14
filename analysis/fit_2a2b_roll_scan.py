import sys, numpy as np
sys.path.insert(0,'analysis'); import ansys_forward as af
from scipy.spatial.transform import Rotation as Rot
from scipy.optimize import minimize
SC='/tmp/claude-0/-home-user-TrinityMEG/8fea0924-981d-5d8c-bcd3-7d068c7590dc/scratchpad/an'
z=np.load('data/J_2a_2b.npz'); xs,ys,zs=z['xs'],z['ys'],z['zs']; J=np.stack([z['Jx'],z['Jy'],z['Jz']],-1).astype(np.float64); mag=np.linalg.norm(J,axis=-1)
nz=np.array(np.nonzero(mag>0)).T
print('2a-2b saline extent x %.0f..%.0f y %.0f..%.0f z %.0f..%.0f mm'%(xs[nz[:,0].min()]*1e3,xs[nz[:,0].max()]*1e3,ys[nz[:,1].min()]*1e3,ys[nz[:,1].max()]*1e3,zs[nz[:,2].min()]*1e3,zs[nz[:,2].max()]*1e3))
c=np.array([0.0704,-0.0354,0.0208]); Q=J.reshape(-1,3).sum(0)*1e-9; qh=Q/np.linalg.norm(Q)
# electrode axis: the segment pair sits at one level; the |J| cloud is elongated along the shaft direction? test: PCA of top voxels
i=np.array(np.unravel_index(np.argsort(mag.ravel())[::-1][:3000],mag.shape)).T; pts=np.stack([xs[i[:,0]],ys[i[:,1]],zs[i[:,2]]],1)
u,s,vt=np.linalg.svd(pts-pts.mean(0),full_matrices=False); print('top-|J| cloud principal axes',np.round(vt,2),'sv',np.round(s,3),' Q dir',np.round(qh,2),' Q.y=%.2f'%qh[1])
# current between contacts: flux of J through the plane perpendicular to Q through c (1 mm slab)
X,Y,Z=np.meshgrid(xs,ys,zs,indexing='ij'); dist=(X-c[0])*qh[0]+(Y-c[1])*qh[1]+(Z-c[2])*qh[2]
slab=np.abs(dist)<0.5e-3; I=np.abs(np.sum(J[slab]@qh)*1e-6)   # A (1 mm^2 per voxel column ~ approx)
print('current through the mid-plane perpendicular to Q: %.2f mA at 1 V  (2a-3a export: 5.97 mA)'%(I*1e3))
# --- forward against the new recording, tip fixed, shaft along Ansys y (as in the 2a-3a model), tip on +y side of level 2
o=np.load(f'{SC}/out_new.npz'); T=o['dev_head_t']; R,t=T[:3,:3],T[:3,3]; m=o['mags']; S=o['pos_dev'][m]@R.T+t; N=o['nrm_dev'][m]@R.T; b=o['pk'][m]; names=list(o['names'][m])
tip=np.array([-0.03403,-0.04753,-0.08755]); shaft=af.R_FIT@np.array([0,1.,0])
G,Qe=af.load_elements('data/J_2a_2b.npz')
def fw(Rm,tip_ans,tp,scale):
    P=(G-tip_ans)@Rm.T+tp; Qr=Qe@Rm.T; out=np.zeros(len(S))
    for k in range(len(S)):
        r=S[k]-P; out[k]=np.sum(np.cross(Qr,r)/np.linalg.norm(r,axis=1)[:,None]**3@N[k])
    return 1e-7*out*scale
scale=7.5e-3/I
sc=lambda bm:(np.corrcoef(bm,b)[0,1],np.dot(bm,b)/np.dot(bm,bm))
for lab,tip_ans in [('tip at +y of level 2',c+np.array([0,0.00275,0])),('tip at -y of level 2',c-np.array([0,0.00275,0]))]:
    res=[]
    for roll in range(0,360,10):
        Rm=(Rot.from_rotvec(np.radians(roll)*shaft)*Rot.from_matrix(af.R_FIT)).as_matrix(); cc,a=sc(fw(Rm,tip_ans,tip,scale)); res.append((cc,a,roll))
    res.sort(key=lambda x:-x[0]); print(f'{lab}: best corr {res[0][0]:.3f} amp ratio {res[0][1]:.2f} at roll {res[0][2]}; next {[(round(r[0],2),r[2]) for r in res[1:4]]}; worst corr {res[-1][0]:.2f}')
    if 'at +y' in lab: best=res[0]
Rm=(Rot.from_rotvec(np.radians(best[2])*shaft)*Rot.from_matrix(af.R_FIT)).as_matrix(); bm=fw(Rm,c+np.array([0,0.00275,0]),tip,scale)
print('model at given tip (x amp ratio): min %.1f max %.1f pT; measured %.1f..%.1f'%(bm.min()*best[1]*1e12,bm.max()*best[1]*1e12,b.min()*1e12,b.max()*1e12))
print('top sensors (meas, model):',[(names[k],f'{b[k]*1e12:+.1f}',f'{bm[k]*best[1]*1e12:+.1f}') for k in np.argsort(-np.abs(b))[:8]])
np.savez(f'{SC}/fit_2a2b.npz',bm=bm*best[1],b=b,roll=best[2],I=I)
