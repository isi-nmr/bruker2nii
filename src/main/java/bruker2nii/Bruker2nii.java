package bruker2nii;

import bruker_plugin_lib.Bruker;
import bruker_plugin_lib.DataBruker;
import com.ericbarnhill.niftijio.Nifti1Header;
import com.ericbarnhill.niftijio.NiftiVolume;
import com.ericbarnhill.niftijio.tools.IndexIterator;
import org.NifTiMRS.DIM_KEYS;
import org.NifTiMRS.NiftiMRS;
import org.NifTiMRS.Nucleus;
import org.NifTiMRS.TAGS;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.lang3.ArrayUtils;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.factory.Nd4j;
import org.nd4j.linalg.indexing.INDArrayIndex;
import org.nd4j.linalg.indexing.NDArrayIndex;
import org.nd4j.linalg.ops.transforms.Transforms;

import java.io.IOException;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.stream.IntStream;


import static org.nd4j.linalg.ops.transforms.Transforms.sqrt;

public class Bruker2nii {
    String pathBruker;
    String pathNii;
    Bruker bruker;
    DataBruker dataBruker;
    NiftiMRS niftiMRS;

    public Bruker2nii(String pathBruker) {
        this.pathBruker = pathBruker;
        bruker = new Bruker();
        bruker.setPath(Paths.get(pathBruker));
        dataBruker = bruker.getData();


    }
    public void convert(String s, char d) {
//        convert(s, d, false, false);
    }
    public boolean convert(String s, DataType type) {
        return convert(s, 't', false, true, type);
    }

    public boolean convert(String pathNii, char d, boolean b, boolean b1, DataType type) {
        this.pathNii = pathNii;
        switch (type) {
            case MRS: {
                boolean squeez = false;
                int brukerLength = bruker.getDims().length;
                if (bruker.getDims()[brukerLength - 1] == 1)
                    squeez = true;
                if (squeez)
                    brukerLength -= 1;

                int[] niftishape = new int[3 + brukerLength];
                int[] shape = Arrays.stream(bruker.getDims()).mapToInt(i -> (int) i).toArray();
                System.arraycopy(new int[]{1, 1, 1}, 0, niftishape, 0, 3);
                if (squeez)
                    System.arraycopy(shape, 0, niftishape, 3, bruker.getDims().length - 1);
                else System.arraycopy(shape, 0, niftishape, 3, bruker.getDims().length);

                niftiMRS = new NiftiMRS(niftishape);

                INDArray PVM_VoxArrSize = bruker.getJcampdx().getMethod().getINDArray("PVM_VoxArrSize");
                try {
                    for (int i = 1; i < 4; i++) { // what is 4?
                        niftiMRS.getNifti().getHeader1().pixdim[i] = PVM_VoxArrSize.getFloat(i - 1);
                    }
                } catch (Exception e) {
                    e.printStackTrace();
                }
                // here pixdim 5 to 7 must be set
                try {
                    niftiMRS.getNifti().getHeader1().pixdim[4] = 1 / bruker.getJcampdx().getMethod().getINDArray("PVM_SpecSWH").getFloat(0);
                } catch (Exception e) {
                    niftiMRS.getNifti().getHeader1().pixdim[4] = bruker.getJcampdx().getMethod().getINDArray("PVM_SpecMatrix").getFloat(0) / bruker.getJcampdx().getMethod().getINDArray("PVM_SpecAcquisitionTime").getFloat(0);
                }
                ArrayList<int[]> idcs;
                if (squeez)
                    idcs = new IndexIterator().iterate(Arrays.copyOfRange(shape, 1, shape.length - 1));
                else idcs = new IndexIterator().iterate(Arrays.copyOfRange(shape, 1, shape.length));

                for (int[] idc : idcs) {
                    int[] idx = new int[niftishape.length];
                    System.arraycopy(idc, 0, idx, 4, idc.length);
                    for (int i = 0; i < shape[0]; i++) {
                        idx[3] = i;
                        idx[0] = 0;
                        int[] dataidx = Arrays.copyOfRange(idx, 3, idx.length);
                        niftiMRS.getNifti().getData().set(idx, dataBruker.real.getFloat(dataidx));
                        idx[0] = 1;
                        niftiMRS.getNifti().getData().set(idx, dataBruker.imag.getFloat(dataidx));
                    }
                }
                String NUCLEUS = bruker.getJcampdx().getAcqp().getString("NUCLEUS");
                if (NUCLEUS.contains("1H"))
                    niftiMRS.getJson().setResonantNucleus(new String[]{Nucleus.N_1H.toString()});
                Float TrnsFreq = bruker.getJcampdx().getAcqp().getFloat("SFO1");
                niftiMRS.getJson().setSpectrometerFrequency(new Double[]{Double.valueOf(TrnsFreq)});

                if (!squeez) { // what is it really!? avg or coil
                    niftiMRS.getJson().setDim_5(DIM_KEYS.DIM_COIL);
                    niftiMRS.getJson().setDim_5_info("coil");
                    niftiMRS.getJson().setDim_5_header(TAGS.EchoTime, Collections.singletonList(bruker.getJcampdx().getTE(10000)));
                    niftiMRS.getJson().setDim_5_header(TAGS.RepetitionTime, Collections.singletonList(bruker.getJcampdx().getTR(10000)));
                }
                INDArray visuCorePosition = bruker.getJcampdx().getMethod().getINDArray("PVM_VoxArrPosition");
                INDArray visuCoreOrientationReShaped = bruker.getJcampdx().getMethod().getINDArray("PVM_VoxArrGradOrient").get(NDArrayIndex.point(0), NDArrayIndex.all());
                INDArray affine_mat = getAffineMat(type, visuCoreOrientationReShaped,visuCorePosition);
                niftiMRS.getNifti().getHeader1().srow_x = affine_mat.getRow(0).toFloatVector();
                niftiMRS.getNifti().getHeader1().srow_y = affine_mat.getRow(1).toFloatVector();
                niftiMRS.getNifti().getHeader1().srow_z = affine_mat.getRow(2).toFloatVector();
                niftiMRS.getNifti().getHeader1().descrip = new StringBuffer("Converted by JBruker2nii API");


                niftiMRS.getNifti().getHeader1().qform_code = 0;
                niftiMRS.getNifti().getHeader1().sform_code = 1;


                try {
                    niftiMRS.write(pathNii, b, b1);
                    return true;
                } catch (IOException e) {
                    e.printStackTrace();
                    return false;
                }
            }
            case MRSI: {
                int[] shape_nifti = new int[3 + 1];
                int[] shape_bruker = Arrays.stream(bruker.getDims()).mapToInt(i -> (int) i).toArray();
                shape_nifti[0] = shape_bruker[1];
                shape_nifti[1] = shape_bruker[2];
                if (shape_bruker.length == 3) {
                    shape_nifti[2] = 1;
                    shape_nifti[3] = shape_bruker[0];
                } else
                    shape_nifti[2] = 1;
                shape_nifti[3] = shape_bruker[0];
                niftiMRS = new NiftiMRS(shape_nifti);
                INDArray PVM_VoxArrSize = bruker.getJcampdx().getMethod().getINDArray("PVM_VoxArrSize");
                INDArray VisuCoreExtent = bruker.getJcampdx().getVisu_pars().getINDArray("VisuCoreExtent");
                INDArray VisuCoreSize = bruker.getJcampdx().getVisu_pars().getINDArray("VisuCoreSize");
                for (int i = 1; i < 3; i++) {
                    niftiMRS.getNifti().getHeader1().pixdim[i] = VisuCoreExtent.getFloat(i) / VisuCoreSize.getFloat(i);
                }
                niftiMRS.getNifti().getHeader1().pixdim[3] = bruker.getJcampdx().getVisu_pars().getFloat("VisuCoreFrameThickness");
                niftiMRS.getNifti().getHeader1().pixdim[4] = VisuCoreExtent.getFloat(0) / VisuCoreSize.getFloat(0);
                ArrayList VisuCoreUnits = bruker.getJcampdx().getVisu_pars().getArrayList("VisuCoreUnits");
//                double[][][] CSI2nii = new double[(int) bruker.getDims()[1]][(int) bruker.getDims()[2]][(int) bruker.getDims()[0]];
//                NiftiHeader header = new Nifti1Header((int) bruker.getDims()[1], (int) bruker.getDims()[1], 1, (int) bruker.getDims()[0]);
                niftiMRS.getNifti().getHeader1().t_unit_code = Nifti1Header.NIFTI_UNITS_SEC;
                if (VisuCoreUnits.contains("mm"))
                    niftiMRS.getNifti().getHeader1().xyz_unit_code = Nifti1Header.NIFTI_UNITS_MM;
                else if (VisuCoreUnits.contains("m"))
                    niftiMRS.getNifti().getHeader1().xyz_unit_code = Nifti1Header.NIFTI_UNITS_METER;

                for (int i = 0; i < 2*bruker.getDims()[1]; i=i+2) {
                    for (int j = 0; j < bruker.getDims()[2]; j++) {
                        INDArrayIndex[] indx = {NDArrayIndex.all(), NDArrayIndex.point(i/2), NDArrayIndex.point(j)};
                        double[] tempR = dataBruker.real.get(indx).toDoubleVector();
                        double[] tempI = dataBruker.imag.get(indx).toDoubleVector();
                        int[] idx = new int[]{i, j, 0, 0};
                        for (int t = 0; t < bruker.getDims()[0]; t++) {
                            idx[3] = t;
                            idx[0] = i;
                            niftiMRS.getNifti().getData().set(idx, tempR[t]);
                            idx[0] = i + 1;
                            niftiMRS.getNifti().getData().set(idx, tempI[t]);
                        }

                    }
                }
                INDArray visuCoreOrientation = bruker.getJcampdx().getVisu_pars().getINDArray("VisuCoreOrientation");
                INDArray visuCoreOrientationReShaped = visuCoreOrientation.get(NDArrayIndex.point(0), NDArrayIndex.all()).reshape(new int[]{3, 3});
                INDArray visuCorePosition = bruker.getJcampdx().getVisu_pars().getINDArray("VisuCorePosition").getRow(0);
                INDArray affine_mat = getAffineMat(type, visuCoreOrientationReShaped,visuCorePosition);
                // if sample upside down
                //            q mat
                INDArray transpos = affine_mat.get(new INDArrayIndex[]{NDArrayIndex.interval(0, 3), NDArrayIndex.point(3)});
                INDArray RZS = affine_mat.get(new INDArrayIndex[]{NDArrayIndex.interval(0, 3), NDArrayIndex.interval(0, 3)});
                INDArray zoom = sqrt(RZS.mul(RZS).sum(0));
                INDArray rotate = RZS.divRowVector(zoom);
                // qfactor
//            int nRows = rotate.rows();
//            int nColumns = rotate.columns();
//            INDArray S = Nd4j.zeros(1, nRows);
//            INDArray U = Nd4j.zeros(nRows, nRows);
//            INDArray V = Nd4j.zeros(nColumns, nColumns);
//            Nd4j.getBlasWrapper().lapack().gesvd(rotate, S, U, V);
//            INDArray PR = U.mmul(V);
//            // error proning
//            PR.reshape(new int[]{1,9});
//          s mat
                niftiMRS.getNifti().getHeader1().srow_x = affine_mat.getRow(0).toFloatVector();
                niftiMRS.getNifti().getHeader1().srow_y = affine_mat.getRow(1).toFloatVector();
                niftiMRS.getNifti().getHeader1().srow_z = affine_mat.getRow(2).toFloatVector();
                String NUCLEUS = bruker.getJcampdx().getAcqp().getString("NUCLEUS");
                if (NUCLEUS.contains("1H"))
                    niftiMRS.getJson().setResonantNucleus(new String[]{Nucleus.N_1H.toString()});
                Float TrnsFreq = bruker.getJcampdx().getAcqp().getFloat("SFO1");
                niftiMRS.getJson().setSpectrometerFrequency(new Double[]{Double.valueOf(TrnsFreq)});
                niftiMRS.getNifti().getHeader1().descrip = new StringBuffer("Converted by JBruker2nii API");
//                    float VisuCoreDataSlope = bruker.getJcampdx().getVisu_pars().getINDArray("VisuCoreDataSlope").getFloat(0);
//                    float VisuCoreDataOffs = bruker.getJcampdx().getVisu_pars().getINDArray("VisuCoreDataOffs").getFloat(0);
//                niftiMRS.getNifti().getHeader1().scl_slope = VisuCoreDataSlope;
//                niftiMRS.getNifti().getHeader1().scl_slope = VisuCoreDataOffs;
                    String slice_order = bruker.getJcampdx().getMethod().getString("PVM_ObjOrderScheme");

                    switch (slice_order) {
                        case "User_defined_slice_scheme" :
                            niftiMRS.getNifti().getHeader1().slice_code = 0;
                            break;
                        case "Sequential" :
                            niftiMRS.getNifti().getHeader1().slice_code = 1;
                            break;
                        case "Reverse_sequential" :
                            niftiMRS.getNifti().getHeader1().slice_code = 2;
                            break;
                        case "Interlaced" :
                            niftiMRS.getNifti().getHeader1().slice_code = 3;
                            break;
                        case "Reverse_interlacesd" :
                            niftiMRS.getNifti().getHeader1().slice_code = 4;
                            break;
                        case "Angiopraphy" :
                            niftiMRS.getNifti().getHeader1().slice_code = 0;
                            break;
                    }
                niftiMRS.getNifti().getHeader1().qform_code = 0;
                niftiMRS.getNifti().getHeader1().sform_code = 1;
                try {
                    niftiMRS.write(pathNii, b, b1);
                    return true;
                } catch (IOException e) {
                    e.printStackTrace();
                    return false;
                }

            }
            case IMAGE: {

                INDArray visuCoreOrientation = bruker.getJcampdx().getVisu_pars().getINDArray("VisuCoreOrientation");
                ArrayList VisuCoreSlicePacksDef = null;
                try {
                    VisuCoreSlicePacksDef = bruker.getJcampdx().getVisu_pars().getArrayList("VisuCoreSlicePacksDef");
                } catch (Exception e) {
                    e.printStackTrace();
                }
                int numOfSlice;
                if (VisuCoreSlicePacksDef == null) {
                    numOfSlice = 1;
                } else {
                    numOfSlice = Integer.valueOf(((Double) VisuCoreSlicePacksDef.get(1)).intValue());
                }
                if (numOfSlice == 1) {
                    int[] shape_nifti = new int[3];
                    int[] shape_bruker = Arrays.stream(bruker.getCplxDims()).mapToInt(i -> (int) i).toArray();
                    shape_nifti[0] = shape_bruker[0];
                    shape_nifti[1] = shape_bruker[1];
                    if (shape_bruker.length>3)
                        dataBruker.real = Nd4j.squeeze(dataBruker.real,-1);
                    if (shape_bruker.length >= 3)
                        shape_nifti[2] = shape_bruker[2];
                    else
                        shape_nifti[2] = 1;
                    int[] dims = new int[]{(int) shape_nifti[0], (int) shape_nifti[1], shape_nifti[2]};
                    Nifti1Header header = new Nifti1Header(dims);
                    INDArray PVM_VoxArrSize = bruker.getJcampdx().getMethod().getINDArray("PVM_VoxArrSize");
                    INDArray VisuCoreExtent = bruker.getJcampdx().getVisu_pars().getINDArray("VisuCoreExtent");
                    INDArray VisuCoreSize = bruker.getJcampdx().getVisu_pars().getINDArray("VisuCoreSize");
                    for (int i = 1; i < 3; i++) {
                        header.pixdim[i] = VisuCoreExtent.getFloat(i - 1) / VisuCoreSize.getFloat(i - 1);
                    }
                    header.pixdim[3] = bruker.getJcampdx().getVisu_pars().getFloat("VisuCoreFrameThickness");
                    ArrayList VisuCoreUnits = bruker.getJcampdx().getVisu_pars().getArrayList("VisuCoreUnits");

                    header.t_unit_code = Nifti1Header.NIFTI_UNITS_SEC;
                    if (VisuCoreUnits.contains("mm"))
                        header.xyz_unit_code = Nifti1Header.NIFTI_UNITS_MM;
                    else if (VisuCoreUnits.contains("m"))
                        header.xyz_unit_code = Nifti1Header.NIFTI_UNITS_METER;
                    NiftiVolume nifti = new NiftiVolume(header);
                    for (int i = 0; i < shape_nifti[0]; i++) {
                        for (int j = 0; j < shape_nifti[1]; j++) {
                            for (int k = 0; k < shape_nifti[2]; k++) {
                                INDArrayIndex[] indx = {NDArrayIndex.point(i), NDArrayIndex.point(j), NDArrayIndex.point(k)};
                                double tempR = dataBruker.real.get(indx).getDouble();
                                int[] idx = new int[]{i, j, k};
                                nifti.getData().set(idx, tempR);
                            }
                        }
                    }

                    INDArray visuCorePosition = bruker.getJcampdx().getVisu_pars().getINDArray("VisuCorePosition").getRow(0);
                    INDArray visuCoreOrientationReShaped = visuCoreOrientation.get(NDArrayIndex.point(0), NDArrayIndex.all()).reshape(new int[]{3, 3});
                    INDArray affine_mat = getAffineMat(type, visuCoreOrientationReShaped,visuCorePosition);
                    nifti.getHeader1().srow_x = affine_mat.getRow(0).toFloatVector();
                    nifti.getHeader1().srow_y = affine_mat.getRow(1).toFloatVector();
                    nifti.getHeader1().srow_z = affine_mat.getRow(2).toFloatVector();
                    nifti.getHeader1().descrip = new StringBuffer("Converted by JBruker2nii API");
                    float VisuCoreDataSlope = bruker.getJcampdx().getVisu_pars().getINDArray("VisuCoreDataSlope").getFloat(0);
                    float VisuCoreDataOffs = bruker.getJcampdx().getVisu_pars().getINDArray("VisuCoreDataOffs").getFloat(0);
                    nifti.getHeader1().scl_slope = VisuCoreDataSlope;
                    nifti.getHeader1().scl_slope = VisuCoreDataOffs;
                    String slice_order = bruker.getJcampdx().getMethod().getString("PVM_ObjOrderScheme");

                    switch (slice_order) {
                        case "User_defined_slice_scheme" :
                            nifti.getHeader1().slice_code = 0;
                            break;
                        case "Sequential" :
                            nifti.getHeader1().slice_code = 1;
                            break;
                        case "Reverse_sequential" :
                            nifti.getHeader1().slice_code = 2;
                            break;
                        case "Interlaced" :
                            nifti.getHeader1().slice_code = 3;
                            break;
                        case "Reverse_interlacesd" :
                            nifti.getHeader1().slice_code = 4;
                            break;
                        case "Angiopraphy" :
                            nifti.getHeader1().slice_code = 0;
                            break;
                    }
                    nifti.getHeader1().qform_code = 0;
                    nifti.getHeader1().sform_code = 1;

                    //                    data_slp = get_value(visu_pars, 'VisuCoreDataSlope')
//                    data_off = get_value(visu_pars, 'VisuCoreDataOffs')
                    try {
                        String extention = FilenameUtils.getExtension(pathNii);
                        if(extention.isEmpty()) {
                            pathNii = pathNii + ".nii";
                        }
                        nifti.write(pathNii);
                        return true;
                    } catch (IOException e) {
                        e.printStackTrace();
                        return false;
                    }
                } else {
                    for (int slice = 0; slice < numOfSlice; slice++) {
                        int[] shape_nifti = new int[3];
                        int[] shape_bruker = Arrays.stream(bruker.getCplxDims()).mapToInt(i -> (int) i).toArray();
                        shape_nifti[0] = shape_bruker[0];
                        shape_nifti[1] = shape_bruker[1];
                        shape_nifti[2] = 1;
                        int[] dims = new int[]{(int) shape_nifti[0], (int) shape_nifti[1], shape_nifti[2]};
                        Nifti1Header header = new Nifti1Header(dims);

                        INDArray PVM_VoxArrSize = bruker.getJcampdx().getMethod().getINDArray("PVM_VoxArrSize");
                        INDArray VisuCoreExtent = bruker.getJcampdx().getVisu_pars().getINDArray("VisuCoreExtent");
                        INDArray VisuCoreSize = bruker.getJcampdx().getVisu_pars().getINDArray("VisuCoreSize");
                        for (int i = 1; i < 3; i++) {
                            header.pixdim[i] = VisuCoreExtent.getFloat(i - 1) / VisuCoreSize.getFloat(i - 1);
                        }
                        header.pixdim[3] = bruker.getJcampdx().getVisu_pars().getFloat("VisuCoreFrameThickness");
                        ArrayList VisuCoreUnits = bruker.getJcampdx().getVisu_pars().getArrayList("VisuCoreUnits");

                        header.t_unit_code = Nifti1Header.NIFTI_UNITS_SEC;
                        if (VisuCoreUnits.contains("mm"))
                            header.xyz_unit_code = Nifti1Header.NIFTI_UNITS_MM;
                        else if (VisuCoreUnits.contains("m"))
                            header.xyz_unit_code = Nifti1Header.NIFTI_UNITS_METER;
                        NiftiVolume nifti = new NiftiVolume(header);
                        for (int i = 0; i < shape_nifti[0]; i++) {
                            for (int j = 0; j < shape_nifti[1]; j++) {
                                INDArrayIndex[] indx = {NDArrayIndex.point(i), NDArrayIndex.point(j), NDArrayIndex.point(slice)};
                                double tempR = dataBruker.real.get(indx).getDouble();
                                int[] idx = new int[]{i, j};
                                nifti.getData().set(idx, tempR);
                            }
                        }

                        INDArray visuCorePosition = bruker.getJcampdx().getVisu_pars().getINDArray("VisuCorePosition").getRow(slice);
                        INDArray visuCoreOrientationReShaped = visuCoreOrientation.get(NDArrayIndex.point(slice), NDArrayIndex.all()).reshape(new int[]{3, 3});
                        INDArray affine_mat = getAffineMat(type, visuCoreOrientationReShaped,visuCorePosition);
                        nifti.getHeader1().srow_x = affine_mat.getRow(0).toFloatVector();
                        nifti.getHeader1().srow_y = affine_mat.getRow(1).toFloatVector();
                        nifti.getHeader1().srow_z = affine_mat.getRow(2).toFloatVector();
                        nifti.getHeader1().descrip = new StringBuffer("Converted by JBruker2nii API");
                        String[] splittedPath = pathNii.split("\\.");
                        String numberedPath = splittedPath[0] + "_" + String.valueOf(slice) + "." + splittedPath[1];
                        float VisuCoreDataSlope = bruker.getJcampdx().getVisu_pars().getINDArray("VisuCoreDataSlope").getFloat(0);
                        float VisuCoreDataOffs = bruker.getJcampdx().getVisu_pars().getINDArray("VisuCoreDataOffs").getFloat(0);
                        nifti.getHeader1().scl_slope = 1;
                        nifti.getHeader1().scl_slope = 0;
                        String slice_order = bruker.getJcampdx().getMethod().getString("PVM_ObjOrderScheme");
                        nifti.getHeader1().datatype = Nifti1Header.NIFTI_TYPE_FLOAT64;
                        switch (slice_order) {
                            case "User_defined_slice_scheme" :
                                nifti.getHeader1().slice_code = 0;
                                break;
                            case "Sequential" :
                                nifti.getHeader1().slice_code = 1;
                                break;
                            case "Reverse_sequential" :
                                nifti.getHeader1().slice_code = 2;
                                break;
                            case "Interlaced" :
                                nifti.getHeader1().slice_code = 3;
                                break;
                            case "Reverse_interlacesd" :
                                nifti.getHeader1().slice_code = 4;
                                break;
                            case "Angiopraphy" :
                                nifti.getHeader1().slice_code = 0;
                                break;
                        }
                        nifti.getHeader1().qform_code = 0;
                        nifti.getHeader1().sform_code = 1;
                        try {
                            nifti.write(numberedPath);
                            return true;
                        } catch (IOException e) {
                            e.printStackTrace();
                            return false;
                        }
                    }

                }
            }
        }

            return false;
        }





//        if (!bruker.isRaw()) {
//            if (isSVS()) {
////                int sigLen = bruker.getJcampdx().getVisu_pars().getINDArray("VisuAcqSize").getInt(0);
//                if (bruker.isIR()) {
//                    int sigLen = bruker.getJcampdx().getVisu_pars().getINDArray("VisuAcqSize").getInt(0);
//                    double[] data2nii = new double[sigLen];
//                    int[] niftishape = new int[3 + bruker.getDims().length];
//                    int[] shape = Arrays.stream(bruker.getDims()).mapToInt(i -> (int) i).toArray();
//                    System.arraycopy(new int[]{1, 1, 1}, 0, niftishape, 0, 3);
//                    System.arraycopy(shape, 0, niftishape, 3, bruker.getDims().length);
//
//                    niftiMRS = new NiftiMRS(niftishape);
//                    INDArray PVM_VoxArrSize = bruker.getJcampdx().getMethod().getINDArray("PVM_VoxArrSize");
//                    for (int i = 1; i < 4; i++) {
//                        niftiMRS.getNifti().header.pixdim[i] = PVM_VoxArrSize.getFloat(i - 1);
//                    }
//                    // here pixdim 5 to 7 must be set
//                    niftiMRS.getNifti().header.pixdim[4] = bruker.getJcampdx().getVisu_pars().getINDArray("VisuCoreExtent").getFloat(0)
//                            / bruker.getJcampdx().getVisu_pars().getINDArray("VisuCoreSize").getFloat(0);
//                    ArrayList<int[]> idcs = new IndexIterator().iterate(Arrays.copyOfRange(shape, 1, shape.length));
//                    for (int[] idc : idcs) {
//                        int[] idx = new int[niftishape.length];
//                        System.arraycopy(idc, 0, idx, 4, idc.length);
//                        for (int i = 0; i < shape[0]; i++) {
//                            idx[3] = i;
//                            idx[0] = 0;
//                            int[] dataidx = Arrays.copyOfRange(idx, 3, idx.length);
//                            niftiMRS.getNifti().data.set(idx, dataBruker.real.getFloat(dataidx));
//                            idx[0] = 1;
//                            niftiMRS.getNifti().data.set(idx, dataBruker.imag.getFloat(dataidx));
//                        }
//                    }
////                    double[][] data2nii = new double[2][sigLen];
////                    data2nii[0] = dataBruker.real.getColumn(0).toDoubleVector();
////                    data2nii[1] = dataBruker.imag.getColumn(0).toDoubleVector();
////                    // check the possibility of repitted measure!
////                    int[] niftishape = new int[]{1, 1, 1, (int) data2nii[0].length};
////                    niftiMRS = new NiftiMRS(niftishape);
////                    INDArray PVM_VoxArrSize = bruker.getJcampdx().getMethod().getINDArray("PVM_VoxArrSize");
////                    for (int i = 1; i < 4; i++) {
////                        niftiMRS.getNifti().header.pixdim[i] = PVM_VoxArrSize.getFloat(i - 1);
////                    }
////                    if (d == 'f') {
////                        niftiMRS.getNifti().header.pixdim[4] = bruker.getJcampdx().getVisu_pars().getINDArray("VisuCoreExtent").getFloat(0);
////                        niftiMRS.getNifti().header.xyz_unit_code = NiftiHeader.NIFTI_UNITS_MM;
////                        niftiMRS.getNifti().header.t_unit_code = NiftiHeader.NIFTI_UNITS_PPM;
//////                        for (int i = 0; i < data2nii[0].length; i++) {
//////                            for (int j = 0; j < 2; j++) {
////////                                volumeW.data.set(j, 0, 0, i, data2nii[j][i]);
//////                            }
//////                        }
////                        int[] shape = Arrays.stream(bruker.getDims()).mapToInt(i -> (int) i).toArray();
////                        ArrayList<int[]> idcs = new IndexIterator().iterate(Arrays.copyOfRange(shape, 1, shape.length));
////
////                        for (int[] idc : idcs) {
////                            int[] idx = new int[niftishape.length];
////                            System.arraycopy(idc, 0, idx, 4, idc.length);
////                            for (int i = 0; i < shape[0]; i++) {
////                                idx[3] = i;
////                                idx[0] = 0;
////                                int[] dataidx = Arrays.copyOfRange(idx, 3, idx.length);
////                                niftiMRS.getNifti().data.set(idx, dataBruker.real.getFloat(dataidx));
////                                idx[0] = 1;
////                                niftiMRS.getNifti().data.set(idx, dataBruker.imag.getFloat(dataidx));
////                            }
////                        }
//                        String NUCLEUS = bruker.getJcampdx().getAcqp().getString("NUCLEUS");
//                        if (NUCLEUS.contains("1H"))
//                            niftiMRS.getJson().ResonantNucleus = new String[]{Nucleus.N_1H.toString()};
//                        Float TrnsFreq = bruker.getJcampdx().getAcqp().getFloat("SFO1");
//                        niftiMRS.getJson().SpectrometerFrequency = new Double[]{Double.valueOf(TrnsFreq)};
////                volumeW.header.scl_slope = bruker.getJcampdx().getVisu_pars().getINDArray("VisuCoreDataSlope").getFloat(0);
////                volumeW.header.scl_inter = bruker.getJcampdx().getVisu_pars().getINDArray("VisuCoreDataOffs").getFloat(0);
//
////                volumeW.header.cal_max = bruker.getJcampdx().getVisu_pars().getFloat("VisuCoreDataMax");
////                volumeW.header.cal_min = bruker.getJcampdx().getVisu_pars().getFloat("VisuCoreDataMin");
//                        niftiMRS.getJson().setDim_5(DIM_KEYS.DIM_MEAS);
//                        niftiMRS.getJson().setDim_5_info("repetition");
//                        niftiMRS.getJson().setDim_5_header(TAGS.EchoTime, Collections.singletonList(bruker.getJcampdx().getTE(10000)));
//                        niftiMRS.getJson().setDim_5_header(TAGS.RepetitionTime, Collections.singletonList(bruker.getJcampdx().getTR(10000)));
//
//                        niftiMRS.getNifti().header.descrip = new StringBuffer("Converted by JBruker2nii API");
//                    }
////                    else if (d == 't') {
////                        FastFourierTransformer fft = new FastFourierTransformer(DftNormalization.STANDARD);
////                        Complex[] complexdata = TransformUtils.createComplexArray(data2nii);
////                        Complex[] timedata = fft.transform(complexdata, TransformType.INVERSE);
////                        data2nii = TransformUtils.createRealImaginaryArray(timedata);
////                        niftiMRS.getNifti().header.xyz_unit_code = NiftiHeader.NIFTI_UNITS_MM;
////                        niftiMRS.getNifti().header.t_unit_code = NiftiHeader.NIFTI_UNITS_SEC;
////                        niftiMRS.getNifti().header.pixdim[4] = (float) (1 / bruker.getJcampdx().getSW(4000));
////                        for (int i = 0; i < data2nii[0].length; i++) {
////                            for (int j = 0; j < 2; j++) {
////                                niftiMRS.getNifti().data.set(j, 0, 0, i, data2nii[j][i]);
////                            }
////                        }
////                    }
//                }
////                else {
////                    double[] data2nii = new double[sigLen];
////                    data2nii = dataBruker.real.getColumn(0).toDoubleVector();
//////                data2nii[1]  = dataBruker.imag.getColumn(0).toDoubleVector();
////                    niftiHeader = new NiftiHeader(1, 1, 1, (int) data2nii.length);
////                    niftiHeader.datatype = NiftiHeader.NIFTI_TYPE_FLOAT64;
////                    volumeW = new NiftiVolume(niftiHeader);
////                    INDArray PVM_VoxArrSize = bruker.getJcampdx().getMethod().getINDArray("PVM_VoxArrSize");
////                    for (int i = 1; i < 4; i++) {
////                        volumeW.header.pixdim[i] = PVM_VoxArrSize.getFloat(i - 1);
////                    }
////                    volumeW.header.pixdim[4] = bruker.getJcampdx().getVisu_pars().getINDArray("VisuCoreExtent").getFloat(0)
////                            / bruker.getJcampdx().getVisu_pars().getINDArray("VisuCoreSize").getFloat(0);
////                    for (int i = 0; i < data2nii.length; i++) {
//////                        volumeW.data.set(0, 0, 0, i, data2nii[i]);
////                    }
////
////                }
////                volumeW.header.scl_slope = bruker.getJcampdx().getVisu_pars().getINDArray("VisuCoreDataSlope").getFloat(0);
////                volumeW.header.scl_inter = bruker.getJcampdx().getVisu_pars().getINDArray("VisuCoreDataOffs").getFloat(0);
////
////                volumeW.header.cal_max = bruker.getJcampdx().getVisu_pars().getFloat("VisuCoreDataMax");
////                volumeW.header.cal_min = bruker.getJcampdx().getVisu_pars().getFloat("VisuCoreDataMin");
////                volumeW.header.descrip = new StringBuffer("Converted by JBruker2nii API");
//// to do affine matrix + 1r+1d + time
//
////            } else {
////
////                // if ( 3d and 4d )
////                double[][][] CSI2nii = new double[(int) bruker.getDims()[1]][(int) bruker.getDims()[2]][(int) bruker.getDims()[0]];
////                NiftiHeader header = new NiftiHeader((int) bruker.getDims()[1], (int) bruker.getDims()[1], 1, (int) bruker.getDims()[0]);
////                header.datatype = NiftiHeader.NIFTI_TYPE_COMPLEX64;
////                volumeW = new NiftiVolume(header);
////
////                INDArray VisuCoreExtent = bruker.getJcampdx().getVisu_pars().getINDArray("VisuCoreExtent");
////                INDArray VisuCoreSize = bruker.getJcampdx().getVisu_pars().getINDArray("VisuCoreSize");
////
////                for (int i = 1; i < 3; i++) {
////                    volumeW.header.pixdim[i] = VisuCoreExtent.getFloat(i) / VisuCoreSize.getFloat(i);
////                }
////                volumeW.header.pixdim[3] = bruker.getJcampdx().getVisu_pars().getFloat("VisuCoreFrameThickness");
////                volumeW.header.pixdim[4] = VisuCoreExtent.getFloat(0) / VisuCoreSize.getFloat(0);
////                ArrayList VisuCoreUnits = bruker.getJcampdx().getVisu_pars().getArrayList("VisuCoreUnits");
////
////                if (VisuCoreUnits.contains("[ppm]"))
////                    volumeW.header.t_unit_code = NiftiHeader.NIFTI_UNITS_PPM;
////                if (VisuCoreUnits.contains("mm"))
////                    volumeW.header.xyz_unit_code = NiftiHeader.NIFTI_UNITS_MM;
////                else if (VisuCoreUnits.contains("m"))
////                    volumeW.header.xyz_unit_code = NiftiHeader.NIFTI_UNITS_METER;
////
////
////                volumeW.header.intent_code = 2001; // each voxel has time_series
////                double[] tempR;
////                double[] tempI;
////                for (int i = 0; i < bruker.getDims()[1]; i++) {
////                    for (int j = 0; j < bruker.getDims()[2]; j=j+1) {
////                        INDArrayIndex[] indx = {NDArrayIndex.all(), NDArrayIndex.point(i), NDArrayIndex.point(j)};
////                        tempR = dataBruker.real.get(indx).toDoubleVector();
////
////                        for (int t = 0; t < bruker.getDims()[0]; t++) {
//////                            volumeW.data.set(j, i, 0, t, tempR[t]);
////                        }
////
////                    }
////                }
//
//        } else if (bruker.isRaw()) {
//            if (isSVS()) {
//                int sigLen = bruker.getJcampdx().getVisu_pars().getINDArray("VisuAcqSize").getInt(0);
//                double[] data2nii = new double[sigLen];
//                int[] niftishape = new int[3 + bruker.getDims().length];
//                int[] shape = Arrays.stream(bruker.getDims()).mapToInt(i -> (int) i).toArray();
//                System.arraycopy(new int[]{1, 1, 1}, 0, niftishape, 0, 3);
//                System.arraycopy(shape, 0, niftishape, 3, bruker.getDims().length);
//
//                niftiMRS = new NiftiMRS(niftishape);
//                INDArray PVM_VoxArrSize = bruker.getJcampdx().getMethod().getINDArray("PVM_VoxArrSize");
//                for (int i = 1; i < 4; i++) {
//                    niftiMRS.getNifti().header.pixdim[i] = PVM_VoxArrSize.getFloat(i - 1);
//                }
//                // here pixdim 5 to 7 must be set
//                niftiMRS.getNifti().header.pixdim[4] = bruker.getJcampdx().getVisu_pars().getINDArray("VisuCoreExtent").getFloat(0)
//                        / bruker.getJcampdx().getVisu_pars().getINDArray("VisuCoreSize").getFloat(0);
//                ArrayList<int[]> idcs = new IndexIterator().iterate(Arrays.copyOfRange(shape, 1, shape.length));
//                for (int[] idc : idcs) {
//                    int[] idx = new int[niftishape.length];
//                    System.arraycopy(idc, 0, idx, 4, idc.length);
//                    for (int i = 0; i < shape[0]; i++) {
//                        idx[3] = i;
//                        idx[0] = 0;
//                        int[] dataidx = Arrays.copyOfRange(idx, 3, idx.length);
//                        niftiMRS.getNifti().data.set(idx, dataBruker.real.getFloat(dataidx));
//                        idx[0] = 1;
//                        niftiMRS.getNifti().data.set(idx, dataBruker.imag.getFloat(dataidx));
//                    }
//                }
//                String NUCLEUS = bruker.getJcampdx().getAcqp().getString("NUCLEUS");
//                if (NUCLEUS.contains("1H"))
//                    niftiMRS.getJson().ResonantNucleus = new String[]{Nucleus.N_1H.toString()};
//                Float TrnsFreq = bruker.getJcampdx().getAcqp().getFloat("SFO1");
//                niftiMRS.getJson().SpectrometerFrequency = new Double[]{Double.valueOf(TrnsFreq)};
////                volumeW.header.scl_slope = bruker.getJcampdx().getVisu_pars().getINDArray("VisuCoreDataSlope").getFloat(0);
////                volumeW.header.scl_inter = bruker.getJcampdx().getVisu_pars().getINDArray("VisuCoreDataOffs").getFloat(0);
//
////                volumeW.header.cal_max = bruker.getJcampdx().getVisu_pars().getFloat("VisuCoreDataMax");
////                volumeW.header.cal_min = bruker.getJcampdx().getVisu_pars().getFloat("VisuCoreDataMin");
//                niftiMRS.getJson().setDim_5(DIM_KEYS.DIM_MEAS);
//                niftiMRS.getJson().setDim_5_info("repetition");
//                niftiMRS.getJson().setDim_5_header(TAGS.EchoTime, Collections.singletonList(bruker.getJcampdx().getTE(10000)));
//                niftiMRS.getJson().setDim_5_header(TAGS.RepetitionTime, Collections.singletonList(bruker.getJcampdx().getTR(10000)));
//
//                niftiMRS.getNifti().header.descrip = new StringBuffer("Converted by JBruker2nii API");
//// to do affine matrix + 1r+1d + time
//            } else {
//                // if ( 3d and 4d )
////                double[][][] CSI2nii = new double[(int) bruker.getDims()[1]][(int) bruker.getDims()[2]][(int) bruker.getDims()[0]];
////                NiftiHeader header = new NiftiHeader((int) bruker.getDims()[1], (int) bruker.getDims()[1], 1, (int) bruker.getDims()[0]);
////                header.datatype = NiftiHeader.NIFTI_TYPE_COMPLEX64;
////                volumeW = new NiftiVolume(header);
////
////                INDArray VisuCoreExtent = bruker.getJcampdx().getVisu_pars().getINDArray("VisuCoreExtent");
////                INDArray VisuCoreSize = bruker.getJcampdx().getVisu_pars().getINDArray("VisuCoreSize");
////
////                for (int i = 1; i < 3; i++) {
////                    volumeW.header.pixdim[i] = VisuCoreExtent.getFloat(i) / VisuCoreSize.getFloat(i);
////                }
////                volumeW.header.pixdim[3] = bruker.getJcampdx().getVisu_pars().getFloat("VisuCoreFrameThickness");
////                volumeW.header.pixdim[4] = VisuCoreExtent.getFloat(0) / VisuCoreSize.getFloat(0);
////                ArrayList VisuCoreUnits = bruker.getJcampdx().getVisu_pars().getArrayList("VisuCoreUnits");
////
////                if (VisuCoreUnits.contains("[ppm]"))
////                    volumeW.header.t_unit_code = NiftiHeader.NIFTI_UNITS_PPM;
////                if (VisuCoreUnits.contains("mm"))
////                    volumeW.header.xyz_unit_code = NiftiHeader.NIFTI_UNITS_MM;
////                else if (VisuCoreUnits.contains("m"))
////                    volumeW.header.xyz_unit_code = NiftiHeader.NIFTI_UNITS_METER;
////
////
////                volumeW.header.intent_code = 2001; // each voxel has time_series
////                double[] tempR;
////                double[] tempI;
////                for (int i = 0; i < bruker.getDims()[1]; i++) {
////                    for (int j = 0; j < 2*bruker.getDims()[2]; j=j+2) {
////                        INDArrayIndex[] indx = {NDArrayIndex.all(), NDArrayIndex.point(i), NDArrayIndex.point(j/2)};
////                        tempR = dataBruker.real.get(indx).toDoubleVector();
////                        tempI = dataBruker.imag.get(indx).toDoubleVector();
////                        for (int t = 0; t < bruker.getDims()[0]; t++) {
//////                            volumeW.data.set(j, i, 0, t, tempR[t]);
//////                            volumeW.data.set(j+1, i, 0, t, tempI[t]);
////                        }
////
////                    }
////                }}
//            }
//        }



    private boolean isSVS() {
        try {
            if (bruker.getJcampdx().getVisu_pars().getArrayList("VisuCoreDimDesc")
                    .equals(new ArrayList() {{add("spectroscopic");}})) {
                return true;
            }
        } catch (Exception e) {
        }
        return false;
    }

    private INDArray getAffineMat(DataType type, INDArray visuCoreOrientationReShaped, INDArray visuCorePosition) {
        // check slice order


        INDArray resol = null;
        INDArray visuCoreSize = null;
        INDArray visuCoreExtent = null;
        String visuSubjectPosition = bruker.getJcampdx().getVisu_pars().getString("VisuSubjectPosition");

        String subjectType = bruker.getJcampdx().getVisu_pars().getString("VisuSubjectType");
        switch (type) {
            case IMAGE:{
                visuCoreSize = bruker.getJcampdx().getVisu_pars().getINDArray("VisuCoreSize");
                visuCoreExtent = bruker.getJcampdx().getVisu_pars().getINDArray("VisuCoreExtent");
                INDArray resolution = visuCoreExtent.div(visuCoreSize);
                INDArray spatialResolution = null;
                if (visuCoreSize.length()<3) {
                    double VisuCoreFrameThickness = bruker.getJcampdx().getVisu_pars().getFloat("VisuCoreFrameThickness");
//        INDArray spatialResolution = Nd4j.append(resolution.get(NDArrayIndex.interval(1, 3)), 1, VisuCoreFrameThickness, -1);
                    spatialResolution = Nd4j.append(resolution, 1, VisuCoreFrameThickness, -1);
                } else {
                    spatialResolution= resolution;
                }
                INDArray pvm_SPackArrGradOrient = bruker.getJcampdx().getMethod().getINDArray("PVM_SPackArrGradOrient");
                INDArray visuCoreUnits = bruker.getJcampdx().getVisu_pars().getINDArray("VisuCoreUnits");
                INDArray gradient_orient = bruker.getJcampdx().getMethod().getINDArray("PVM_SPackArrGradOrient");
                // add frame count to visuCoreSize (just for image)
                Float VisuCoreDataSlope = bruker.getJcampdx().getVisu_pars().getFloat("VisuCoreDataSlope");
                Float VisuCoreDataOffs = bruker.getJcampdx().getVisu_pars().getFloat("VisuCoreDataOffs");
                // slope Correction
                // offset Correction
                String[] subPos = new String[]{"Head_Prone", "Head_Supine"};
                if (Arrays.stream(subPos).anyMatch(e -> e.contains(visuSubjectPosition))) {
                    System.out.println("Known case ('Head_Prone' or  'Head_Supine' for the parameter 'visu_pars.VisuSubjectPosition.)");
                }

                INDArray axis_orient = Nd4j.zeros(3);
                axis_orient.put(0, Transforms.abs(visuCoreOrientationReShaped).getRow(0).argMax());
                axis_orient.put(1, Transforms.abs(visuCoreOrientationReShaped).getRow(1).argMax());
                axis_orient.put(2, Transforms.abs(visuCoreOrientationReShaped).getRow(2).argMax());
                int subj_orient = axis_orient.getInt(2);
                if (subj_orient == 0 || subj_orient == 2 ) {
                    resol = Nd4j.diag(spatialResolution);
                } else {
                    resol = Nd4j.diag(spatialResolution.mul(Nd4j.createFromArray(new float[]{1,1,-1})));
                }
                break;
            }
            case MRSI: {
                visuCoreSize = bruker.getJcampdx().getVisu_pars().getINDArray("VisuCoreSize").get(NDArrayIndex.interval(1, 3));
                visuCoreExtent = bruker.getJcampdx().getVisu_pars().getINDArray("VisuCoreExtent").get(NDArrayIndex.interval(1, 3));
                INDArray resolution = visuCoreExtent.div(visuCoreSize);
                double VisuCoreFrameThickness = bruker.getJcampdx().getVisu_pars().getFloat("VisuCoreFrameThickness");
//        INDArray spatialResolution = Nd4j.append(resolution.get(NDArrayIndex.interval(1, 3)), 1, VisuCoreFrameThickness, -1);
                INDArray spatialResolution = Nd4j.append(resolution, 1, VisuCoreFrameThickness, -1);
                INDArray pvm_SPackArrGradOrient = bruker.getJcampdx().getMethod().getINDArray("PVM_SPackArrGradOrient");
                INDArray visuCoreUnits = bruker.getJcampdx().getVisu_pars().getINDArray("VisuCoreUnits");
                INDArray gradient_orient = bruker.getJcampdx().getMethod().getINDArray("PVM_SPackArrGradOrient");
                // add frame count to visuCoreSize (just for image)
                Float VisuCoreDataSlope = bruker.getJcampdx().getVisu_pars().getFloat("VisuCoreDataSlope");
                Float VisuCoreDataOffs = bruker.getJcampdx().getVisu_pars().getFloat("VisuCoreDataOffs");
                // slope Correction
                // offset Correction
                String[] subPos = new String[]{"Head_Prone", "Head_Supine"};
                if (Arrays.stream(subPos).anyMatch(e -> e.contains(visuSubjectPosition))) {
                    System.out.println("Known case ('Head_Prone' or  'Head_Supine' for the parameter 'visu_pars.VisuSubjectPosition.)");
                }

                INDArray axis_orient = Nd4j.zeros(3);
                axis_orient.put(0, Transforms.abs(visuCoreOrientationReShaped).getRow(0).argMax());
                axis_orient.put(1, Transforms.abs(visuCoreOrientationReShaped).getRow(1).argMax());
                axis_orient.put(2, Transforms.abs(visuCoreOrientationReShaped).getRow(2).argMax());
                int subj_orient = axis_orient.getInt(2);
                if (subj_orient == 0 || subj_orient == 2 ) {
                    resol = Nd4j.diag(spatialResolution);
                } else {
                    resol = Nd4j.diag(spatialResolution.mul(Nd4j.createFromArray(new float[]{1,1,-1})));
                }
                break;
            }
            case MRS: {
                INDArray spatialResolution = bruker.getJcampdx().getMethod().getINDArray("PVM_VoxArrSize");
                String[] subPos = new String[]{"Head_Prone", "Head_Supine"};
                if (Arrays.stream(subPos).anyMatch(e -> e.contains(visuSubjectPosition))) {
                    System.out.println("Known case ('Head_Prone' or  'Head_Supine' for the parameter 'visu_pars.VisuSubjectPosition.)");
                }
                INDArray axis_orient = Nd4j.zeros(3);
                axis_orient.put(0, Transforms.abs(visuCoreOrientationReShaped).getRow(0).argMax());
                axis_orient.put(1, Transforms.abs(visuCoreOrientationReShaped).getRow(1).argMax());
                axis_orient.put(2, Transforms.abs(visuCoreOrientationReShaped).getRow(2).argMax());
                int subj_orient = axis_orient.getInt(2);
                if (subj_orient == 0 || subj_orient == 2 ) {
                    resol = Nd4j.diag(spatialResolution);
                } else {
                    resol = Nd4j.diag(spatialResolution.mul(Nd4j.createFromArray(new float[]{1,1,-1})));
                }

            }
        }



        INDArray visuCoreOrientationReShaped_trans = visuCoreOrientationReShaped.transpose().mmul(resol);

        INDArray result = Nd4j.eye(4);
        result.put(new INDArrayIndex[]{NDArrayIndex.interval(0, 3), NDArrayIndex.interval(0, 3)}, visuCoreOrientationReShaped_trans);
        result.put(new INDArrayIndex[]{NDArrayIndex.interval(0, 3), NDArrayIndex.point(3)}, visuCorePosition);
//                  Invert the orientation matrix, according to nifti convention and Bruker manual.
        INDArray result_rotated = null;
        switch (visuSubjectPosition) {
            case "Head_Supine":{
                result_rotated = apply_rotate(result, new double[]{0d, 0d, Math.PI});
                break;
            }
            case "Head_Prone":{
                result_rotated = result;
                break;
            }
            case "Head_Left":{
                result_rotated = apply_rotate(result, new double[]{0d, 0d, Math.PI/2});
                break;
            }
            case "Head_Right":{
                result_rotated = apply_rotate(result, new double[]{0d, 0d, -Math.PI/2});
                break;
            }
            case "Foot_Supine":
            case "Tail_Supine":{
                result_rotated = apply_rotate(result, new double[]{Math.PI, 0d,0d });
                break;
            }
            case "Foot_Prone":
            case "Tail_Prone":{
                result_rotated = apply_rotate(result, new double[]{0d, Math.PI,0d});
                break;
            }
            case "Foot_Left":
            case "Tail_Left":{
                result_rotated = apply_rotate(result, new double[]{0d, 0d, Math.PI/2});
                break;
            }
            case "Foot_Right":
            case "Tail_Right":{
                result_rotated = apply_rotate(result, new double[]{0d, 0d, -Math.PI/2});
                break;
            }
            default:
                System.out.println("ERROR: cannot apply visuSubjectPosition");
        }

        try {
            if(!subjectType.equals("Biped")){
                result_rotated = apply_rotate(result_rotated, new double[]{-Math.PI/2, Math.PI, 0});
            }
        } catch (Exception e) {
            System.out.println("subjectType is not defined");
        }


//        InvertMatrix.invert(result, true);
//        // calculate determinant
//        float[][] defMat = new float[][]{{1, 0, 0}, {0, 0, 1}, {0, 1, 0}};
//        INDArray result_orientation =
//                result.get(new INDArrayIndex[]{NDArrayIndex.interval(0, 3), NDArrayIndex.interval(0, 3)}).mmul(Nd4j.createFromArray(defMat));
////        from SAR to ASL
//        if (pivot(result_orientation.getColumn(0).toDoubleVector()) > 0) {
//            result_orientation.getColumn(0).muli(-1);
//        }
//        if (pivot(result_orientation.getColumn(1).toDoubleVector()) > 0) {
//            result_orientation.getColumn(1).muli(-1);
//        }
//        if (pivot(result_orientation.getColumn(2).toDoubleVector()) < 0) {
//            result_orientation.getColumn(2).muli(-1);
//        }
//        result_orientation.mmuli(Nd4j.diag(spatialResolution));
//        result.put(new INDArrayIndex[]{NDArrayIndex.interval(0, 3), NDArrayIndex.interval(0, 3)}, result_orientation);
        // consider_subject_position:
        // keep same det
        return result_rotated;
    }

    public static Double pivot(double[] arr) {
        Double[] doubleArray = ArrayUtils.toObject(arr);
        List<Double> arrList = Arrays.asList(doubleArray);
        int indx = IntStream.range(0, arr.length).reduce((a, b) -> ((double) Math.abs(arrList.get(a))) < ((double) Math.abs(arrList.get(b))) ? b : a).getAsInt();
        return arr[indx];
    }

    public static INDArray calc_eulerangle(INDArray mati) {
        INDArray mat = mati.get(NDArrayIndex.point(0), NDArrayIndex.all(), NDArrayIndex.all());
        double sy = Math.sqrt(mat.getDouble(0, 0) * mat.getDouble(0, 0) + mat.getDouble(1, 0) * mat.getDouble(1, 0));
        double x, y, z;
        if (sy < 1e-6) {
            x = Math.atan2(-mat.getDouble(1, 2), mat.getDouble(1, 1));
            y = Math.atan2(-mat.getDouble(2, 0), sy);
            z = 0;
        } else {
            x = Math.atan2(mat.getDouble(2, 1), mat.getDouble(2, 2));
            y = Math.atan2(-mat.getDouble(2, 0), sy);
            z = Math.atan2(mat.getDouble(1, 0), mat.getDouble(0, 0));
        }
        return Nd4j.create(new double[]{x, y, z});
    }
    public INDArray apply_rotate(INDArray mat, double[] rad){
        INDArray x = Nd4j.createFromArray(new float[]{1, 0, 0, 0, (float) Math.cos(rad[0]), (float) -Math.sin(rad[0]), 0, (float) Math.sin(rad[0]), (float) Math.cos(rad[0])}).reshape( new int[]{3, 3});
        INDArray y = Nd4j.createFromArray(new float[] {(float) Math.cos(rad[1]), 0, (float) Math.sin(rad[1]),0, 1, 0, (float) -Math.sin(rad[1]),0, (float) Math.cos(rad[1])}).reshape( new int[]{3, 3});
        INDArray z = Nd4j.createFromArray(new float[] {(float) Math.cos(rad[2]), (float) -Math.sin(rad[2]), 0, (float) Math.sin(rad[2]), (float) Math.cos(rad[2]), 0, 0,0,1}).reshape( new int[]{3, 3});
        INDArray rotated_mat = z.mmul(y.mmul(x.mmul(mat.get(new INDArrayIndex[]{NDArrayIndex.interval(0, 3), NDArrayIndex.interval(0, 3)}))));
        INDArray rotated_vec = z.mmul(y.mmul(x.mmul(mat.get(new INDArrayIndex[]{NDArrayIndex.interval(0, 3),NDArrayIndex.point(3)}))));
        INDArray result = Nd4j.eye(4);
        result.put(new INDArrayIndex[]{NDArrayIndex.interval(0, 3), NDArrayIndex.interval(0, 3)}, rotated_mat);
        result.put(new INDArrayIndex[]{NDArrayIndex.interval(0, 3), NDArrayIndex.point(3)}, rotated_vec);
        return result;
    }

}
