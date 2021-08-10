//import com.ericbarnhill.niftijio.NiftiHeader;
//import com.ericbarnhill.niftijio.NiftiVolume;
//
//import java.io.IOException;
//import java.io.PrintWriter;
//
//public class Exmple_amir {
//    public static void main(String[] args) throws IOException {
////        NiftiHeader header = new NiftiHeader(1, 1, 1, 10);
//        NiftiHeader header = new NiftiHeader(new int[] {1,1,1,10});
//        header.datatype = NiftiHeader.NIFTI_TYPE_COMPLEX64;
//        NiftiVolume volumeW = new NiftiVolume(header);
//        for (int i = 0; i < 10; i++) {
//            for (int j = 0; j <2 ; j++) {
//                volumeW.data.set(j,0,0,i,i);
//            }
//        }
//
//        volumeW.header.pixdim[4] = (float) 0.00025;
//        volumeW.write("example.nii.gz");
//
//        NiftiVolume volume = NiftiVolume.read("example.nii.gz");
//
//        int nx = volume.header.dim[1];
//        int ny = volume.header.dim[2];
//        int nz = volume.header.dim[3];
//        int dim = volume.header.dim[4];
//
//        if (dim == 0)
//            dim = 1;
//
//        if (args.length == 1)
//        {
//            System.out.println("dims: " + nx + " " + ny + " " + nz + " " + dim);
//            System.out.println("datatype: " + NiftiHeader.decodeDatatype(volume.header.datatype));
//        }
//        else if (args[1].endsWith("txt"))
//        {
//            PrintWriter out = new PrintWriter(args[1]);
//
//            out.println("volume ");
//            out.println("dimensions:");
//            out.println(nx + " " + ny + " " + nz + " " + dim);
//            out.println("data:");
//            for (int d = 0; d < dim; d++)
//                for (int k = 0; k < nz; k++)
//                    for (int j = 0; j < ny; j++)
//                        for (int i = 0; i < nx; i++)
//                            out.println(volume.data.get(i,j,k,d));
//
//            out.println();
//            out.close();
//        }
//        else
//        {
//            volume.write(args[1]);
//        }
//    }
//}
