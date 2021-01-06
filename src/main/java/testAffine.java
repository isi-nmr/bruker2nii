public class testAffine {
    public static void main(String[] args) {
        Bruker2nii bruker2nii = new Bruker2nii("D:\\my-project\\bruker2nifti\\test_data\\bru_banana\\1\\pdata\\1\\2dseq");
        bruker2nii.convert("testAffine.nii.gz");
    }
}
