package bruker2nii;

public class testAffine {
    public static void main(String[] args) {
//        Bruker2nii bruker2nii = new Bruker2nii("D:\\bruker APIs\\brukerRaw\\20190724_114946_BRKRAW_1_1\\5\\pdata\\1\\2dseq");
//        bruker2nii.convert("testAffine5.nii", DataType.IMAGE);
//        Bruker2nii bruker2nii2 = new Bruker2nii("D:\\DATA SETs\\20210118_135716_180121_2mixedPhantoms_Amir_CSI_PRESS_180121_1_1\\T2_TurboRARE(E11)\\pdata\\1\\2dseq");
//        bruker2nii2.convert("rare.nii", DataType.IMAGE);
//        Bruker2nii bruker2nii3 = new Bruker2nii("D:\\DATA SETs\\20210121_085814_21012021_1mixedPhantom_CSI_PRESS_21012021_1_1_1\\CSI_notBad(E14)_30\\fid");
//        bruker2nii3.convert("csi", DataType.MRSI);
//        Bruker2nii bruker2nii4 = new Bruker2nii("D:\\DATA SETs\\SSY-059 ZS, JS, IP, AS\\data\\phantom_IP_AS_Amirs_concentrations\\20200724_125056_MIX4_c2_1_1\\PRESS_1H_hermite(E4)\\fid");
//        bruker2nii4.convert("press", DataType.MRS);
//        Bruker2nii bruker2nii5 = new Bruker2nii("D:\\DATA SETs\\20210809_143251_phantomMIX1_phantomMIX1_1_1\\EPI(E15)\\pdata\\1\\2dseq");
//        bruker2nii5.convert("epi3d.nii", DataType.IMAGE);

//        Bruker2nii test = new Bruker2nii("D:\\DATA SETs\\20210809_143251_phantomMIX1_phantomMIX1_1_1\\T2_TurboRARE(E7)\\pdata\\1\\2dseq");
//        test.convert("test.nii", DataType.IMAGE);
//        Bruker2nii test = new Bruker2nii("D:\\DATA SETs\\20210809_143251_phantomMIX1_phantomMIX1_1_1\\T2_TurboRARE(E9)\\pdata\\1\\2dseq");
//        test.convert("test2.nii", DataType.IMAGE);
//        Bruker2nii test = new Bruker2nii("D:\\DATA SETs\\20210809_143251_phantomMIX1_phantomMIX1_1_1\\CSI(E12)\\fid");
//        test.convert("test4.nii", DataType.MRSI);
        Bruker2nii test = new Bruker2nii("D:\\DATA SETs\\20210809_143251_phantomMIX1_phantomMIX1_1_1\\pavlovaPRESS(E3)\\fid");
        test.convert("epi.nii", DataType.MRS);

    }
}
