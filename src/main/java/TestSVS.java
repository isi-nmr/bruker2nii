import java.io.IOException;

public class TestSVS {
    public static void main(String[] args) throws IOException {
        Bruker2nii bruker2nii = new Bruker2nii("D:\\test fid-a\\FID-A\\exampleData\\Bruker\\sample01_press\\press\\pdata\\1\\1r");
//        Bruker2nii bruker2nii = new Bruker2nii("D:\\DATA SETs\\for test Jbruker\\20200612_094625_lego_phantom_3_1_2\\35\\fid");
//        Bruker2nii bruker2nii = new Bruker2nii("D:\\DATA SETs\\SSY-059 ZS, JS, IP, AS\\data\\phantom_IP_AS\\20200520_151057_phantomNAA_pH_7_2_1_2\\5\\fid");

        bruker2nii.convert("testSVS2", 'f', true,false);

//        Bruker bruker = new Bruker();
//        String pathTemp = "D:\\DATA SETs\\for test Jbruker\\20200612_094625_lego_phantom_3_1_2\\35\\pdata\\1\\2dseq";
//        bruker.setPath(Paths.get(pathTemp));
//        DataBruker data = bruker.getData();
//        double[] realdata0 = data.real.getColumn(0).toDoubleVector();
//        double[] realdata1 = data.real.getColumn(0).toDoubleVector();
//        double[] imagdata0 = data.imag.getColumn(0).toDoubleVector();
//        double[] imagdata1 = data.imag.getColumn(0).toDoubleVector();
//        double[][] data2nii = new double[2][2048];
//        data2nii[0] = realdata0;
//        data2nii[1] = imagdata0;
//
//        NiftiHeader header = new NiftiHeader(data2nii.length/2, 1, 1, (int) data2nii[0].length);
//        header.datatype = NiftiHeader.NIFTI_TYPE_COMPLEX64;
//        NiftiVolume volumeW = new NiftiVolume(header);
//        volumeW.header.pixdim[4] = (float) (1/bruker.getJcampdx().getSW(4000));
//        volumeW.header.xyzt_units = NiftiHeader.NIFTI_UNITS_SEC | NiftiHeader.NIFTI_UNITS_METER;
//        volumeW.header.qform_code = NiftiHeader.NIFTI_XFORM_ALIGNED_ANAT;
//        float quatern_b = 1;
//        float quatern_c = 1;
//        float quatern_d = 1;
//        volumeW.header.quatern = new float[]{quatern_b, quatern_c, quatern_d};
//        volumeW.header.qfac = 1;
//        volumeW.header.intent_code = 2001; // each voxel has time_series
//
//
//        for (int i = 0; i < data2nii[0].length; i++) {
//            for (int j = 0; j < data2nii.length ; j++) {
//                volumeW.data.set(j,0,0,i, data2nii[j][i]);
//            }
//        }
//
//        try {
//            volumeW.write("testSVS.nii.gz");
//        } catch (IOException e) {
//            e.printStackTrace();
//        }

    }
}
