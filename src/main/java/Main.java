import bruker_plugin_lib.Bruker;
import bruker_plugin_lib.DataBruker;
import com.ericbarnhill.niftijio.NiftiHeader;
import com.ericbarnhill.niftijio.NiftiVolume;
import org.nd4j.linalg.indexing.INDArrayIndex;
import org.nd4j.linalg.indexing.NDArrayIndex;

import java.io.IOException;
import java.nio.file.Paths;

public class Main {
    public static void main(String[] args) {
        Bruker bruker = new Bruker();
        String pathTemp = "D:\\DATA SETs\\for test Jbruker\\20200612_094625_lego_phantom_3_1_2\\35\\fid";
        bruker.setPath(Paths.get(pathTemp));
        DataBruker data = bruker.getData();
        double[] realdata0 = data.real.getColumn(0).toDoubleVector();
        double[] realdata1 = data.real.getColumn(0).toDoubleVector();
        double[] imagdata0 = data.imag.getColumn(0).toDoubleVector();
        double[] imagdata1 = data.imag.getColumn(0).toDoubleVector();
        double[][] data2nii = new double[2][2048];
        data2nii[0] = realdata0;
        data2nii[1] = imagdata0;
//        data2nii[2] = realdata1;
//        data2nii[3] = imagdata1;
        NiftiHeader header = new NiftiHeader(data2nii.length/2, 1, 1, (int) data2nii[0].length);
        header.datatype = NiftiHeader.NIFTI_TYPE_COMPLEX64;
        NiftiVolume volumeW = new NiftiVolume(header);
        for (int i = 0; i < data2nii[0].length; i++) {
            for (int j = 0; j < data2nii.length ; j++) {
//                volumeW.data.set(j,0,0,i, data2nii[j][i]);
            }
        }
        volumeW.header.pixdim[4] = (float) (1/bruker.getJcampdx().getSW(4000));
        try {
            volumeW.write("example_1DNMR.nii.gz");
        } catch (IOException e) {
            e.printStackTrace();
        }

        Bruker brukerCSI = new Bruker();
        String pathTempcsi = "D:\\DATA SETs\\for test Jbruker\\20200612_094625_lego_phantom_3_1_2\\30\\fid";
        brukerCSI.setPath(Paths.get(pathTempcsi));
        DataBruker dataCSI = brukerCSI.getData();
        double[][][][] CSI2nii = new double[2][16][16][2048];

        NiftiHeader headerCSI = new NiftiHeader(16, 16, 1, (int) 2048);
        header.datatype = NiftiHeader.NIFTI_TYPE_COMPLEX64;
        NiftiVolume volumeWCSI = new NiftiVolume(header);


        for (int i = 0; i < 16; i++) {
            for (int j = 0; j < 16 ; j++) {
                INDArrayIndex[] indx = { NDArrayIndex.all() , NDArrayIndex.point(i), NDArrayIndex.point(j)};
                CSI2nii[0][i][j] = dataCSI.real.get(indx).toDoubleVector();
                CSI2nii[1][i][j] = dataCSI.imag.get(indx).toDoubleVector();
                for (int t = 0; t < 2048 ; t++) {
//                    volumeW.data.set(j,i,0,t, CSI2nii[0][j][i][t]);
//                    volumeW.data.set(j,i,0,t, CSI2nii[1][j][i][t]);
                }

            }
        }
        volumeW.header.pixdim[4] = (float) (1/brukerCSI.getJcampdx().getSW(4000));
        try {
            volumeW.write("example_CSI.nii.gz");
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
