#! /usr/bin/env python3

# Two-shell phase-reversed diffusion processing using MRTRIX

import re
import os
import pdb
import sys
import glob
import zipfile
import subprocess
import numpy as np

debug_speed_mode = False

os.environ["PATH"] += os.pathsep + '/mrtrix3/bin'
os.environ["PATH"] += os.pathsep + '/usr/local/fsl5/'
os.environ["PATH"] += os.pathsep + '/usr/local/fsl/bin/'
os.environ['FSLDIR'] = '/usr/local/fsl'
# os.environ['FSLDIR'] = '/usr/share/fsl/5.0'
os.environ['FSLOUTPUTTYPE'] = 'NIFTI_GZ'
os.environ['ANTSPATH'] = '/usr/lib/ants'

os.environ['LD_LIBRARY_PATH'] = '/usr/local/fsl/lib'
os.environ['PWD'] = '/flywheel/v0'
os.environ['FLYWHEEL'] = '/flywheel/v0'
os.environ['SHLVL'] = '1'

os.environ["PATH"] += os.pathsep + '/usr/lib/ants'
os.environ["PATH"] += os.pathsep + '/usr/bin'
os.environ["PATH"] += os.pathsep + '/usr/local/bin'
os.environ["PATH"] += os.pathsep + '/usr/local/sbin'
os.environ["PATH"] += os.pathsep + '/usr/sbin'
os.environ["PATH"] += os.pathsep + '/sbin'
os.environ["PATH"] += os.pathsep + '/bin'

# print(os.environ)

work_dir = os.path.join('.', 'work')
output_dir = os.path.join('.', 'output')
# QC_dir = os.path.join(work_dir, 'QC')

os.system('ls -la')
os.makedirs(work_dir, exist_ok=True)
os.system(f'chmod 777 {work_dir}')
os.system('ls -la')

# os.makedirs('work/scratch', exist_ok=True)
# os.makedirs(QC_dir, exist_ok=True)

with open('./output/stdout.txt', 'w') as stdout, open('./output/stderr.txt', 'w') as stderr:
    dwi_zipped_dicoms = [ff.strip() for ff in os.popen(f'find . -type f -maxdepth 3 -regex ".*.dicom.*"').readlines()]
    dwi_zipped_dicoms = [x for x in dwi_zipped_dicoms if os.path.basename(x)[0] != '.']

    [print(x) for x in dwi_zipped_dicoms]

    # UNCOMMENT start
    for dwi_zipped_dicom in dwi_zipped_dicoms:
        shell_dir = os.path.dirname(dwi_zipped_dicom)

        with zipfile.ZipFile(dwi_zipped_dicom, 'r') as zip_ref:
            zip_ref.extractall(shell_dir)

        os.popen(f'ls {shell_dir}')

        output = subprocess.check_output(['dcm2niix', '-w', '1', shell_dir], universal_newlines=True)

        os.popen(f'ls {shell_dir}')
    # UNCOMMENT end

    phen_d = {'AP': 'j', 'RL': 'i'}

    # Prep rest of phase encoding stuff
    _phen_d = dict([(kk[::-1], phen_d[kk] + '-') for kk in phen_d.keys()])
    phen_d_ = phen_d.copy()
    phen_d_.update(_phen_d)

    zeropad_dirs = {0: 'A', 1: 'L', 2: 'I'}

    b0_idx = {}


    def os_find(srch, s_dir='.'):
        result = [ll.strip() for ll in os.popen(f'find {s_dir} -name {srch}').readlines()]
        return result


    def bv_proc(bvalf):
        srcdir = os.path.dirname(bvalf)
        fname = os.path.basename(bvalf)

        if fname[:5] == 'adj_':
            return

        bv = np.loadtxt(bvalf)
        bv[bv < 100] = 0
        np.savetxt(f'{srcdir}/adj_{fname}', bv.T, newline='     ', fmt='%i')
        return bv.max()


    # Find DWIs
    dwi_niftis = [ff.strip() for ff in os.popen(f'find . -maxdepth 3 -name *DWI*nii').readlines()]
    phase_enc_assign = {}  # Dict for PE assignments
    bmax_assign = {}  # Dict for Bmax assignments
    imps = []
    shell_dirs = {}

    # Figure out PE for each series
    for dwi_nifti in dwi_niftis:
        sfn = re.split('[_.]', dwi_nifti)  # Split the file name
        sift_pe = list(set(phen_d_.keys()).intersection(sfn))
        # Sift the filename to grab the PE direction
        if len(sift_pe) == 0:
            continue
        else:
            pe = sift_pe[0]

        ser = sfn[-2]  # Grab the series number

        try:
            s = int(ser)
        except:
            ser = sfn[-3]

        phase_enc_assign[ser] = pe

        shell_dir = os.path.dirname(dwi_nifti)
        shell_dirs[ser] = shell_dir

    # Prep b-vals and get b-max for each series
    for ss in phase_enc_assign.keys():
        # os.system("dcm2niix -w 1 %s" % dcmdirs[ss])
        bvalfs = os_find("*bval", s_dir=shell_dirs[ss])
        for bv in bvalfs:
            if os.path.split(bv)[-1].startswith('adj'):
                os.system('rm ' + bv)

        bvalfs = os_find("*bval", s_dir=shell_dirs[ss])
        bvalf = [bb for bb in bvalfs if os.path.basename(bb)[:5] != 'adj_'][0]  ###
        adj_bvalf = bvalf.replace(shell_dirs[ss] + '/', shell_dirs[ss] + '/' + 'adj_')
        max_b = int(bv_proc(bvalf))
        bmax_assign[ss] = max_b
        b0_idx[ss] = np.where(np.genfromtxt(adj_bvalf) < 100)[0][0]

    # Do converts,  denoises, and list
    for ss in phase_enc_assign.keys():
        # impmif = "dwi_%s_%s.mif" % (bmax_ass[ss],pe_ass[ss])
        bvecf = os_find("*.bvec", s_dir=shell_dirs[ss])[0]
        abvalf = os_find("adj_*bval", s_dir=shell_dirs[ss])[0]
        # Check the grid spacing
        nifti = os.path.splitext(bvecf)[0] + '.nii'
        # Analyze dimenions
        dims = os.popen(f'mrinfo -size {nifti}').readlines()[0].strip().split()[:3]
        need_rs = False
        zp_opts = []

        for i_, d_ in enumerate(dims):
            if int(d_) % 2 != 0:
                need_rs = True
                zp_opts.append(f'-{zeropad_dirs[i_]} 1')

        if need_rs:
            os.system(f'3dZeropad -overwrite -prefix {nifti} {" ".join(zp_opts)} {nifti} ')
        # Convert NIFTI to MIF

        mif = os.path.join(shell_dirs[ss], f'dwi_{bmax_assign[ss]}_{phase_enc_assign[ss]}_{ss}.mif')

        # UNCOMMENT start
        os.system(f'mrconvert {nifti} {mif} -fslgrad {bvecf} {abvalf} '
                  f'-set_property PhaseEncodingDirection {phen_d_[phase_enc_assign[ss]]} '
                  f'-set_property EchoTime 0.098 '
                  f'-set_property TotalReadoutTime 0.0764 '
                  f'-force')
        # UNCOMMENT end

        ##List commands for processing

    proc_cmds = []

    proc_cmds.append('set -x')
    proc_cmds.append('set -e')

    # timing
    proc_cmds.append('start=$SECONDS')
    # for denoises
    for ss in phase_enc_assign.keys():
        impmif = os.path.join(shell_dirs[ss], f'dwi_{bmax_assign[ss]}_{phase_enc_assign[ss]}_{ss}.mif')
        impmif_out = impmif.replace("dwi_", "dn_dwi_")
        proc_cmds.append(f'dwidenoise {impmif} {impmif_out} -force')
        imps.append(impmif_out)
    # timing
    proc_cmds.append('time=$((SECONDS-start)); echo dwidenoise time: $time sec.')

    # Do cat and preproc
    catname = os.path.join('.', 'output',
                           f'dwi_{list(phase_enc_assign.keys())[0]}_{max(bmax_assign.values())}+{len(bmax_assign.values())}sh.mif')

    if len(phase_enc_assign) > 1:
        proc_cmds.append(f'mrcat {" ".join(imps)} {catname} -axis 3 -force')
    else:
        os.system(f'ln -s {imps[0]} {catname}')

    # Extract the b=0 volumes for TOPUP correction
    b0_vols = []
    for ss in phase_enc_assign.keys():
        impmif = os.path.join(shell_dirs[ss], f'dn_dwi_{bmax_assign[ss]}_{phase_enc_assign[ss]}_{ss}.mif')
        impmif_out = impmif.replace('dn_', 'b0_dn_')
        # extract b=0 volume
        idx = b0_idx[ss]
        proc_cmds.append(f'mrconvert -coord 3 {idx} {impmif} {impmif_out} -force')
        b0_vols.append(impmif_out)

    # test if the diffusion scans are properly reversely phase encoded; if yes,
    # use the -se_epi option
    se_epi_str = ''
    dirs = list(set(phase_enc_assign.values()))
    if len(dirs) == 2 and dirs[0] == dirs[1][::-1]:
        se_epi_str = f'-se_epi {catname.replace("dwi", "b0_dn_dwi")} -align_seepi'

    # concatenate the b=0 volumes
    if len(phase_enc_assign) > 1:
        proc_cmds.append(f'mrcat {" ".join(b0_vols)} {catname.replace("dwi", "b0_dn_dwi")} -axis 3 -force')
    else:
        os.system(f'ln -s {b0_vols[0]} {catname.replace("dwi", "b0_dn_dwi")} ')

    # timing
    proc_cmds.append('start=$SECONDS')

    dwipreproc_output = catname.replace('dwi_', 'pp_dwi_')

    eddy_speed_params = '--repol --very_verbose=True'
    if debug_speed_mode:
        eddy_speed_params = '--niter=1 --flm=movement --nvoxhp=100 --repol=False --mbs_niter=1 --dont_peas=False --log_timings=True --very_verbose=True'

    proc_cmds.append(f'echo Entering dwifslpreproc')

    proc_cmds.append(f'tree {work_dir};')
    proc_cmds.append(f'tree {output_dir};')

    proc_cmds.append(
        f'dwifslpreproc {catname} {dwipreproc_output} -rpe_header -force {se_epi_str} -nocleanup -scratch {work_dir} -debug -topup_options " --verbose " -eddy_options " {eddy_speed_params} " 1> ./output/stdout.txt 2> ./output/stderr.txt')
    # -eddyqc_all {output_dir}
    # timing
    proc_cmds.append('time=$((SECONDS-start)); echo dwifslpreproc time: $time sec.')

    proc_cmds.append(f'echo dwifslpreproc DONE')

    proc_cmds.append(f'tree {work_dir};')
    proc_cmds.append(f'tree {output_dir};')

    dwipreproc_output_nii = dwipreproc_output.replace(".mif", ".nii")
    proc_cmds.append(f'mrconvert {dwipreproc_output} {dwipreproc_output_nii} -force')

    fslmask_output = os.path.join('.', 'output', 'fslmask')
    proc_cmds.append(f'bet2 {dwipreproc_output_nii} {fslmask_output} -m -f 0.05 -v')

    bias_dwipreproc_output_nii = dwipreproc_output_nii.replace("pp_", "bpp_")
    proc_cmds.append(
        f'dwibiascorrect fsl {dwipreproc_output} {bias_dwipreproc_output_nii}.gz -mask {fslmask_output}_mask.nii.gz -force')

    # # rename files
    # proc_cmds.append(f'find ./ -iname \'bvals*\' -exec mv \'{{}}\' {output_dir}/ \;')
    # proc_cmds.append(f'find ./ -iname \'bvecs*\' -exec mv \'{{}}\' {output_dir}/ \;')
    # # proc_cmds.append(f'find {work_dir}/ -iname \'{os.path.basename(bias_dwipreproc_output_nii)}.gz\' -exec mv \'{{}}\' {output_dir}/dwi_2sh.nii.gz \;')
    # # proc_cmds.append(f'find {work_dir}/ -iname \'fslmask_mask.nii.gz\' -exec mv \'{{}}\' {output_dir}/b0_dwi_brain_mask_2sh.nii.gz \;')

    # # Remove large intermediate files
    # proc_cmds.append(f'rm {output_dir}/*.mif')
    # proc_cmds.append(f'rm {output_dir}/*.nii')

    proc_cmds.append(f'echo Before bvals and bvecs are copied')

    proc_cmds.append(f'tree {work_dir};')
    proc_cmds.append(f'tree {output_dir};')

    proc_cmds.append(f'find {work_dir} -iname \'bvals*\';')
    proc_cmds.append(f'find {work_dir} -iname \'bvecs*\';')

    proc_cmds.append(f'mv `find {work_dir} -iname \'bvals*\'`  {output_dir};')
    proc_cmds.append(f'mv `find {work_dir} -iname \'bvecs*\'`  {output_dir};')

    proc_cmds.append(f'echo After bvals and bvecs are copied')

    proc_cmds.append(f'tree {work_dir};')
    proc_cmds.append(f'tree {output_dir};')

    # proc_cmds.append(f'tar -zcvf {output_dir}/QC.tar.gz ./QC')
    # proc_cmds.append(f'tar -zcvf {output_dir}/scratch.tar.gz {work_dir}/scratch')

    proc_cmd_out = '\n'.join(proc_cmds)
    proc_file_name = os.path.join('.', 'output', 'proc_this.dti.sh')
    open(proc_file_name, 'w').write(proc_cmd_out)

    os.system(f'chmod +x {proc_file_name}')
    os.system(proc_file_name)

    # bet2 pp_dwi.nii fslmask -m -f 0.05
    # dwibiascorrect pp_dwi.mif bpp_dwi.nii -mask fslmask_mask.nii.gz -fsl -force

    a = 1