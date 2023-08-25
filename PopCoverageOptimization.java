import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Scanner;

public class PopCoverageOptimization {

    public static void main(String[] args)
    {
//         String data = "HLA-A*11:01\tG13D\tVVVGAGDVGK\n" +
//                 "HLA-A*68:01\tG13D\tVVVGAGDVGK\n" +
//                 "HLA-A*03:01\tG13D\tVVVGAGDVGK\n" +
//                 "HLA-A*02:03\tG12D\tKLVVVGADGV\n" +
//                 "HLA-A*02:03\tG13D\tKLVVVGAGDV\n" +
//                 "HLA-A*30:01\tG13D\tVVVGAGDVGK\n" +
//                 "HLA-A*02:01\tG12D\tKLVVVGADGV\n" +
//                 "HLA-A*02:01\tG13D\tKLVVVGAGDV\n" +
//                 "HLA-A*02:06\tG12D\tKLVVVGADGV\n" +
//                 "HLA-A*02:06\tG13C\tKLVVVGAGC\n" +
//                 "HLA-A*02:03\tG13C\tKLVVVGAGC\n" +
//                 "HLA-A*02:01\tG13C\tKLVVVGAGC\n" +
//                 "HLA-A*02:06\tG13D\tKLVVVGAGDV\n" +
//                 "HLA-A*31:01\tG13D\tVVVGAGDVGK\n" +
//                 "HLA-A*02:03\tG13S\tKLVVVGAGS\n" +
//                 "HLA-B*15:01\tG12S\tKLVVVGASG\n" +
//                 "HLA-A*02:01\tG13S\tKLVVVGAGS\n" +
//                 "HLA-A*31:01\tG12D+G13R\tKLVVVGADR\n" +
//                 "HLA-A*68:01\tG12D+G13C\tVVVGADCVGK\n" +
//                 "HLA-A*68:01\tG12C+G13D\tVVVGACDVGK\n" +
//                 "HLA-A*68:01\tG12D+G13D\tVVVGADDVGK\n" +
//                 "HLA-A*11:01\tG12D+G13C\tVVVGADCVGK\n" +
//                 "HLA-A*11:01\tG12C+G13D\tVVVGACDVGK\n" +
//                 "HLA-A*11:01\tG12D+G13D\tVVGADDVGK\n" +
//                 "HLA-A*11:01\tG12D+G13D\tVVVGADDVGK\n" +
//                 "HLA-A*02:06\tG12D+G13S\tLVVVGADSV\n" +
//                 "HLA-A*02:03\tG12D+G13S\tKLVVVGADSV\n" +
//                 "HLA-A*02:03\tG12C+G13S\tKLVVVGACSV\n" +
//                 "HLA-A*68:02\tG12D+G13S\tLVVVGADSV\n" +
//                 "HLA-A*03:01\tG12D+G13C\tVVVGADCVGK\n" +
//                 "HLA-A*03:01\tG12C+G13D\tVVVGACDVGK\n" +
//                 "HLA-A*68:01\tG12D+G13D\tVVGADDVGK\n" +
//                 "HLA-A*68:02\tG12D+G13R\tLVVVGADRV\n" +
//                 "HLA-A*02:06\tG12D+G13R\tLVVVGADRV\n" +
//                 "HLA-A*03:01\tG12D+G13D\tVVVGADDVGK\n" +
//                 "HLA-A*03:01\tG12D+G13D\tVVGADDVGK\n" +
//                 "HLA-A*03:01\tG12D+G13R\tKLVVVGADR\n" +
//                 "HLA-A*02:01\tG12D+G13S\tKLVVVGADSV\n" +
//                 "HLA-A*02:03\tG12D+G13S\tLVVVGADSV\n" +
//                 "HLA-A*02:01\tG12C+G13S\tKLVVVGACSV\n" +
//                 "HLA-A*02:06\tG12D+G13S\tKLVVVGADSV\n" +
//                 "HLA-A*68:02\tG12D+G13C\tLVVVGADCV\n" +
//                 "HLA-A*02:06\tG12D+G13C\tLVVVGADCV\n" +
//                 "HLA-A*02:03\tG12D+G13C\tKLVVVGADCV\n" +
//                 "HLA-A*02:06\tG12C+G13S\tKLVVVGACSV\n" +
//                 "HLA-A*68:01\tG12D+G13R\tKLVVVGADR\n" +
//                 "HLA-A*02:03\tG12D+G13R\tLVVVGADRV\n" +
//                 "HLA-B*51:01\tG12D+G13S\tLVVVGADSV\n" +
//                 "HLA-A*02:01\tG12D+G13S\tLVVVGADSV\n" +
//                 "HLA-A*33:01\tG12D+G13R\tKLVVVGADR\n" +
//                 "HLA-A*30:01\tG12D+G13D\tVVGADDVGK\n" +
//                 "HLA-A*30:01\tG12D+G13C\tVVVGADCVGK\n" +
//                 "HLA-A*30:01\tG12C+G13D\tVVVGACDVGK\n" +
//                 "HLA-A*02:01\tG12D+G13C\tKLVVVGADCV\n" +
//                 "HLA-A*02:03\tG12D+G13D\tKLVVVGADDV\n" +
//                 "HLA-A*02:01\tG12D+G13R\tLVVVGADRV\n" +
//                 "HLA-B*51:01\tG12D+G13R\tLVVVGADRV\n" +
//                 "HLA-A*11:01\tG12D+G13R\tKLVVVGADR\n" +
//                 "HLA-A*31:01\tG12C+G13D\tVVVGACDVGK\n" +
//                 "HLA-A*31:01\tG12D+G13C\tVVVGADCVGK\n" +
//                 "HLA-A*02:03\tG12D+G13C\tLVVVGADCV\n" +
//                 "HLA-A*02:01\tG12D+G13C\tLVVVGADCV\n" +
//                 "HLA-A*02:01\tG12D+G13D\tKLVVVGADDV\n" +
//                 "HLA-A*02:01\tG12D+G13C\tKLVVVGADC\n" +
//                 "HLA-A*26:01\tG12D+G13S\tLVVVGADSV\n" +
//                 "HLA-B*35:01\tG12D+G13S\tLVVVGADSV\n" +
//                 "HLA-A*02:03\tG12S+G13D\tKLVVVGASDV\n" +
//                 "HLA-A*02:01\tG12S+G13D\tKLVVVGASDV\n" +
//                 "HLA-A*02:06\tG12S+G13C\tKLVVVGASC\n" +
//                 "HLA-A*02:03\tG12S+G13C\tKLVVVGASC\n" +
//                 "HLA-A*02:01\tG12S+G13C\tKLVVVGASC\n" +
//                 "HLA-A*02:06\tG12S+G13D\tKLVVVGASDV\n" +
//                 "HLA-A*02:03\tG12S+G13S\tKLVVVGASS\n" +
//                 "HLA-A*02:06\tG12S+G13S\tKLVVVGASS\n" +
//                 "HLA-A*02:01\tG12S+G13S\tKLVVVGASS\n" +
//                 "HLA-B*15:01\tG12S+G13C\tKLVVVGASC\n" +
//                 "HLA-B*15:01\tG12S+G13S\tKLVVVGASS\n" +
//                 "HLA-A*32:01\tG12S+G13C\tKLVVVGASC\n" +
//                 "HLA-A*01:01\tE62G\tILDTAGQGEY\n" +
//                 "HLA-A*01:01\tA59T+Q61R+E62G\tILDTTGRGEY\n" +
//                 "HLA-A*30:02\tA59T+Q61R+E62G\tILDTTGRGEY\n" +
//                 "HLA-A*30:02\tE62G\tILDTAGQGEY\n" +
//                 "HLA-B*15:01\tE62G\tILDTAGQGEY\n" +
//                 "HLA-B*15:01\tA59T+Q61R+E62G\tILDTTGRGEY\n" +
//                 "HLA-A*68:02\tQ61L\tDTAGLEEYSA\n" +
//                 "HLA-A*01:01\tE62G\tLDTAGQGEY\n" +
//                 "HLA-A*01:01\tA59T+Q61H+E62G\tLDTTGHGEY\n" +
//                 "HLA-A*02:01\tQ61L\tLLDILDTAGL\n" +
//                 "HLA-B*35:01\tE62G\tILDTAGQGEY\n" +
//                 "HLA-A*68:01\tQ61R\tLDILDTAGR\n" +
//                 "HLA-B*35:01\tA59T+Q61L+E62G\tTGLGEYSAM\n" +
//                 "HLA-A*30:02\tA59T+Q61H+E62G\tLDTTGHGEY\n" +
//                 "HLA-A*68:02\tA59T+Q61R+E62G\tDTTGRGEYSA\n" +
//                 "HLA-A*30:02\tE62G\tLDTAGQGEY\n" +
//                 "HLA-A*02:03\tQ61L\tLLDILDTAGL\n" +
//                 "HLA-B*35:01\tA59T+Q61R+E62G\tILDTTGRGEY\n" +
//                 "HLA-A*02:06\tQ61L\tLLDILDTAGL\n" +
//                 "HLA-B*51:01\tA59T+Q61L+E62G\tTGLGEYSAM\n" +
//                 "HLA-A*68:01\tQ61R\tLLDILDTAGR\n" +
//                 "HLA-A*26:01\tE62G\tILDTAGQGEY\n" +
//                 "HLA-B*40:01\tQ61L\tLDILDTAGL\n" +
//                 "HLA-A*03:01\tA59T+Q61R+E62G\tILDTTGRGEY\n" +
//                 "HLA-A*33:01\tQ61R\tLDILDTAGR\n" +
//                 "HLA-B*08:01\tA59T+Q61L+E62G\tTGLGEYSAM\n" +
//                 "HLA-A*01:01\tQ61H\tLLDILDTAGH\n" +
//                 "HLA-A*02:06\tA59T+Q61H+E62G\tILDTTGHGE\n" +
//                 "HLA-A*03:01\tE62G\tILDTAGQGEY\n" +
//                 "HLA-B*44:02\tA59T+Q61H+E62G\tLDTTGHGEY\n" +
//                 "HLA-A*26:01\tA59T+Q61R+E62G\tILDTTGRGEY\n" +
//                 "HLA-A*33:01\tQ61R\tLLDILDTAGR\n" +
//                 "HLA-B*44:03\tA59T+Q61H+E62G\tLDTTGHGEY\n" +
//                 "HLA-A*02:01\tA59T+Q61H+E62G\tILDTTGHGE\n" +
//                 "HLA-A*31:01\tQ61R\tLLDILDTAGR\n" +
//                 "HLA-B*44:03\tE62G\tLDTAGQGEY\n" +
//                 "HLA-B*44:02\tE62G\tLDTAGQGEY\n" +
//                 "HLA-B*57:01\tA59T+Q61R+E62G\tILDTTGRGEY\n" +
//                 "HLA-A*01:01\tA59T+Q61H+E62G\tILDTTGHGE\n" +
//                 "HLA-B*15:01\tA59T+Q61L+E62G\tTGLGEYSAM\n" +
//                 "HLA-B*15:01\tA59T+Q61H+E62G\tLDTTGHGEY\n" +
//                 "HLA-A*26:01\tA59T+Q61H+E62G\tLDTTGHGEY\n" +
//                 "HLA-A*02:06\tE62G\tILDTAGQGEY\n" +
//                 "HLA-B*35:01\tA59T+Q61H+E62G\tLDTTGHGEY\n" +
//                 "HLA-B*15:01\tE62G\tLDTAGQGEY\n" +
//                 "HLA-B*57:01\tE62G\tILDTAGQGEY\n" +
//                 "HLA-B*07:02\tA59T+Q61L+E62G\tTGLGEYSAM\n" +
//                 "HLA-A*02:01\tE62G\tILDTAGQGEY\n" +
//                 "HLA-A*26:01\tE62G\tLDTAGQGEY\n" +
//                 "HLA-B*58:01\tA59T+Q61R+E62G\tILDTTGRGEY\n" +
//                 "HLA-B*35:01\tE62G\tLDTAGQGEY\n" +
//                 "HLA-A*30:02\tA59T+Q61L+E62G\tTGLGEYSAM\n" +
//                 "HLA-A*01:01\tE62G\tILDTAGQGE\n" +
//                 "HLA-B*58:01\tE62G\tILDTAGQGEY\n" +
//                 "HLA-A*03:01\tQ61R\tLLDILDTAGR\n" +
//                 "HLA-A*11:01\tE62G\tILDTAGQGEY\n" +
//                 "HLA-A*01:01\tQ61L\tLLDILDTAGL\n" +
//                 "HLA-B*53:01\tE62G\tILDTAGQGEY\n" +
//                 "HLA-B*53:01\tA59T+Q61L+E62G\tTGLGEYSAM\n" +
//                 "HLA-A*31:01\tQ61R\tLDILDTAGR\n" +
//                 "HLA-A*02:01\tA59T+Q61R+E62G\tILDTTGRGEY\n" +
//                 "HLA-A*01:01\tQ61R\tLLDILDTAGR\n" +
//                 "HLA-A*32:01\tA59T+Q61R+E62G\tILDTTGRGEY\n" +
//                 "HLA-A*68:02\tA59T+Q61L+E62G\tTGLGEYSAM\n" +
//                 "HLA-A*02:01\tA59T+Q61R+E62G\tILDTTGRGE\n" +
//                 "HLA-B*53:01\tA59T+Q61R+E62G\tILDTTGRGEY\n" +
//                 "HLA-B*44:02\tA59T+Q61R+E62G\tILDTTGRGEY\n" +
//                 "HLA-A*02:01\tE62G\tILDTAGQGE\n" +
//                 "HLA-A*11:01\tA59T+Q61R+E62G\tILDTTGRGEY\n" +
//                 "HLA-B*44:02\tE62G\tILDTAGQGEY\n" +
//                 "HLA-A*32:01\tE62G\tILDTAGQGEY\n" +
//                 "HLA-A*26:01\tA59T+Q61L+E62G\tTGLGEYSAM\n" +
//                 "HLA-B*44:03\tE62G\tILDTAGQGEY\n" +
//                 "HLA-A*11:01\tQ61R\tLLDILDTAGR\n" +
//                 "HLA-A*26:01\tQ61L\tLDILDTAGL\n" +
//                 "HLA-A*26:01\tQ61L\tDTAGLEEYSA\n" +
//                 "HLA-A*26:01\tA59T\tDILDTTGQE\n" +
//                 "HLA-B*53:01\tA59T+Q61H+E62G\tLDTTGHGEY\n" +
//                 "HLA-B*44:03\tA59T+Q61R+E62G\tILDTTGRGEY\n" +
//                 "HLA-B*44:02\tQ61L\tLDILDTAGL\n" +
//                 "HLA-A*01:01\tA59T+Q61H\tILDTTGHEEY\n" +
//                 "HLA-A*01:01\tA59T+E62G\tILDTTGQGEY\n" +
//                 "HLA-A*01:01\tQ61R+E62G\tILDTAGRGEY\n" +
//                 "HLA-A*30:02\tA59T+Q61H\tILDTTGHEEY\n" +
//                 "HLA-A*30:02\tA59T+E62G\tILDTTGQGEY\n" +
//                 "HLA-A*30:02\tQ61R+E62G\tILDTAGRGEY\n" +
//                 "HLA-B*15:01\tA59T+Q61H\tILDTTGHEEY\n" +
//                 "HLA-A*68:02\tQ61L+E62G\tDTAGLGEYSA\n" +
//                 "HLA-B*15:01\tA59T+E62G\tILDTTGQGEY\n" +
//                 "HLA-B*15:01\tQ61R+E62G\tILDTAGRGEY\n" +
//                 "HLA-B*08:01\tA59T+Q61R\tTGREEYSAM\n" +
//                 "HLA-A*68:02\tQ61R+E62G\tDTAGRGEYSA\n" +
//                 "HLA-B*35:01\tA59T+Q61H\tTGHEEYSAM\n" +
//                 "HLA-A*26:01\tA59T+Q61H\tTTGHEEYSAM\n" +
//                 "HLA-B*07:02\tA59T+Q61R\tTGREEYSAM\n" +
//                 "HLA-A*01:01\tQ61H+E62G\tLDTAGHGEY\n" +
//                 "HLA-A*01:01\tQ61L+E62G\tLDTAGLGEY\n" +
//                 "HLA-B*08:01\tA59T+Q61H\tTGHEEYSAM\n" +
//                 "HLA-B*35:01\tA59T+Q61L\tTGLEEYSAM\n" +
//                 "HLA-B*35:01\tA59T+Q61H\tILDTTGHEEY\n" +
//                 "HLA-A*01:01\tQ61R+E62G\tLDTAGRGEY\n" +
//                 "HLA-B*35:01\tA59T+E62G\tILDTTGQGEY\n" +
//                 "HLA-A*68:02\tA59T+Q61L\tDTTGLEEYSA\n" +
//                 "HLA-B*35:01\tA59T+Q61R\tTGREEYSAM\n" +
//                 "HLA-A*68:01\tA59T+Q61R\tLDILDTTGR\n" +
//                 "HLA-B*51:01\tA59T+Q61H\tTGHEEYSAM\n" +
//                 "HLA-B*08:01\tA59T+Q61L\tTGLEEYSAM\n" +
//                 "HLA-A*68:02\tA59T+Q61H\tTTGHEEYSA\n" +
//                 "HLA-A*30:02\tQ61H+E62G\tLDTAGHGEY\n" +
//                 "HLA-B*51:01\tA59T+Q61L\tTGLEEYSAM\n" +
//                 "HLA-A*68:02\tA59T+Q61H\tDTTGHEEYSA\n" +
//                 "HLA-A*30:02\tQ61L+E62G\tLDTAGLGEY\n" +
//                 "HLA-A*30:02\tQ61R+E62G\tLDTAGRGEY\n" +
//                 "HLA-A*26:01\tA59T+Q61L\tTTGLEEYSAM\n" +
//                 "HLA-A*68:02\tA59T+Q61H\tTTGHEEYSAM\n" +
//                 "HLA-B*15:01\tA59T+Q61R\tTGREEYSAM\n" +
//                 "HLA-B*40:01\tA59T+Q61L\tLDILDTTGL\n" +
//                 "HLA-A*68:01\tA59T+Q61R\tLLDILDTTGR\n" +
//                 "HLA-A*68:02\tA59T+Q61L\tTTGLEEYSA\n" +
//                 "HLA-A*03:01\tA59T+Q61H\tILDTTGHEEY\n" +
//                 "HLA-B*51:01\tA59T+Q61R\tTGREEYSAM\n" +
//                 "HLA-B*35:01\tQ61R+E62G\tILDTAGRGEY\n" +
//                 "HLA-B*15:01\tQ61L+E62G\tAGLGEYSAM\n" +
//                 "HLA-A*26:01\tA59T+E62G\tILDTTGQGEY\n" +
//                 "HLA-B*15:01\tA59T+Q61H\tTGHEEYSAM\n" +
//                 "HLA-A*02:06\tA59T+Q61H\tILDTTGHEEY\n" +
//                 "HLA-A*68:02\tA59T+Q61L\tDTTGLEEYS\n" +
//                 "HLA-A*26:01\tA59T+Q61H\tILDTTGHEEY\n" +
//                 "HLA-A*30:01\tA59T+Q61R\tTGREEYSAM\n" +
//                 "HLA-A*03:01\tA59T+E62G\tILDTTGQGEY\n" +
//                 "HLA-A*02:01\tA59T+Q61H\tILDTTGHEEY\n" +
//                 "HLA-A*68:02\tA59T+Q61L\tTTGLEEYSAM\n" +
//                 "HLA-B*07:02\tA59T+Q61H\tTGHEEYSAM\n" +
//                 "HLA-B*58:01\tA59T+Q61H\tILDTTGHEEY\n" +
//                 "HLA-B*57:01\tA59T+Q61H\tILDTTGHEEY\n" +
//                 "HLA-A*26:01\tA59T+Q61H\tTGHEEYSAM\n" +
//                 "HLA-A*03:01\tQ61R+E62G\tILDTAGRGEY\n" +
//                 "HLA-A*31:01\tA59T+Q61R\tLLDILDTTGR\n" +
//                 "HLA-B*44:02\tQ61H+E62G\tLDTAGHGEY\n" +
//                 "HLA-B*44:03\tQ61H+E62G\tLDTAGHGEY\n" +
//                 "HLA-A*33:01\tA59T+Q61R\tLLDILDTTGR\n" +
//                 "HLA-A*33:01\tA59T+Q61R\tLDILDTTGR\n" +
//                 "HLA-B*53:01\tA59T+Q61H\tILDTTGHEEY\n" +
//                 "HLA-A*68:02\tA59T+Q61H\tDTTGHEEYS\n" +
//                 "HLA-A*26:01\tQ61L+E62G\tLDTAGLGEY\n" +
//                 "HLA-A*26:01\tQ61R+E62G\tILDTAGRGEY\n" +
//                 "HLA-A*30:02\tA59T+Q61R\tTGREEYSAM\n" +
//                 "HLA-A*30:02\tQ61L+E62G\tAGLGEYSAM\n" +
//                 "HLA-B*35:01\tQ61L+E62G\tAGLGEYSAM\n" +
//                 "HLA-B*44:03\tQ61R+E62G\tLDTAGRGEY\n" +
//                 "HLA-A*02:06\tQ61H+E62G\tILDTAGHGE\n" +
//                 "HLA-A*02:01\tA59T+E62G\tILDTTGQGEY\n" +
//                 "HLA-A*32:01\tA59T+Q61H\tILDTTGHEEY\n" +
//                 "HLA-A*02:06\tA59T+E62G\tILDTTGQGEY\n" +
//                 "HLA-B*15:01\tQ61L+E62G\tLDTAGLGEY\n" +
//                 "HLA-A*03:01\tA59T+Q61R\tLLDILDTTGR\n" +
//                 "HLA-B*44:02\tQ61R+E62G\tLDTAGRGEY\n" +
//                 "HLA-B*08:01\tQ61L+E62G\tAGLGEYSAM\n" +
//                 "HLA-A*26:01\tA59T+Q61R\tTGREEYSAM\n" +
//                 "HLA-A*26:01\tQ61H+E62G\tLDTAGHGEY\n" +
//                 "HLA-B*15:01\tQ61H+E62G\tLDTAGHGEY\n" +
//                 "HLA-A*11:01\tA59T+Q61H\tILDTTGHEEY\n" +
//                 "HLA-B*57:01\tQ61R+E62G\tILDTAGRGEY\n" +
//                 "HLA-A*68:01\tQ61L+E62G\tDTAGLGEYSA\n" +
//                 "HLA-A*02:06\tQ61L+E62G\tAGLGEYSAM\n" +
//                 "HLA-A*30:02\tA59T+Q61H\tTGHEEYSAM\n" +
//                 "HLA-A*02:01\tQ61H+E62G\tILDTAGHGE\n" +
//                 "HLA-B*15:01\tA59T+Q61L\tTGLEEYSAM\n" +
//                 "HLA-B*35:01\tQ61H+E62G\tLDTAGHGEY\n" +
//                 "HLA-B*15:01\tQ61R+E62G\tLDTAGRGEY\n" +
//                 "HLA-B*57:01\tA59T+Q61H\tTTGHEEYSAM\n" +
//                 "HLA-B*57:01\tA59T+Q61L\tTTGLEEYSAM\n" +
//                 "HLA-A*26:01\tQ61L+E62G\tDTAGLGEYSA\n" +
//                 "HLA-B*07:02\tQ61L+E62G\tAGLGEYSAM\n" +
//                 "HLA-B*57:01\tA59T+E62G\tILDTTGQGEY\n" +
//                 "HLA-A*01:01\tQ61H+E62G\tILDTAGHGE\n" +
//                 "HLA-B*15:01\tA59T+Q61H\tTTGHEEYSAM\n" +
//                 "HLA-B*53:01\tA59T+Q61L\tTGLEEYSAM\n" +
//                 "HLA-A*26:01\tQ61R+E62G\tLDTAGRGEY\n" +
//                 "HLA-B*44:02\tQ61L+E62G\tLDTAGLGEY\n" +
//                 "HLA-A*26:01\tQ61R+E62G\tDTAGRGEYSA\n" +
//                 "HLA-B*35:01\tQ61R+E62G\tLDTAGRGEY\n" +
//                 "HLA-B*51:01\tQ61L+E62G\tAGLGEYSAM\n" +
//                 "HLA-B*44:03\tQ61L+E62G\tLDTAGLGEY\n" +
//                 "HLA-A*02:06\tA59T+Q61L\tTGLEEYSAM\n" +
//                 "HLA-A*01:01\tA59T+E62G\tILDTTGQGE\n" +
//                 "HLA-A*01:01\tA59T+Q61R\tLLDILDTTGR\n" +
//                 "HLA-B*53:01\tA59T+E62G\tILDTTGQGEY\n" +
//                 "HLA-A*11:01\tA59T+E62G\tILDTTGQGEY\n" +
//                 "HLA-B*58:01\tA59T+E62G\tILDTTGQGEY\n" +
//                 "HLA-B*53:01\tA59T+Q61H\tTGHEEYSAM\n" +
//                 "HLA-B*07:02\tA59T+Q61L\tTGLEEYSAM\n" +
//                 "HLA-B*58:01\tQ61R+E62G\tILDTAGRGEY\n" +
//                 "HLA-A*68:02\tA59T+Q61H\tTGHEEYSAM\n" +
//                 "HLA-B*35:01\tQ61L+E62G\tLDTAGLGEY\n" +
//                 "HLA-B*44:02\tA59T+Q61H\tILDTTGHEEY\n" +
//                 "HLA-B*51:01\tQ61R+E62G\tTAGRGEYSA\n" +
//                 "HLA-A*68:02\tA59T+Q61L\tTGLEEYSAM\n" +
//                 "HLA-B*35:01\tA59T+Q61H\tTTGHEEYSAM\n" +
//                 "HLA-A*68:01\tQ61R+E62G\tDTAGRGEYSA\n" +
//                 "HLA-A*31:01\tA59T+Q61R\tLDILDTTGR\n" +
//                 "HLA-B*08:01\tA59T+Q61H\tILDTTGHEEY\n" +
//                 "HLA-B*44:03\tA59T+Q61H\tILDTTGHEEY\n" +
//                 "HLA-A*68:02\tQ61R+E62G\tTAGRGEYSA\n" +
//                 "HLA-A*01:01\tA59T+Q61H\tTTGHEEYSAM\n" +
//                 "HLA-A*01:01\tA59T+Q61H\tTTGHEEYSA\n" +
//                 "HLA-B*35:01\tA59T+Q61L\tTTGLEEYSAM\n" +
//                 "HLA-A*01:01\tA59T+Q61L\tTTGLEEYSA\n" +
//                 "HLA-A*02:01\tA59T+E62G\tILDTTGQGE\n" +
//                 "HLA-A*01:01\tA59T+Q61L\tTTGLEEYSAM\n" +
//                 "HLA-A*26:01\tA59T+Q61L\tTGLEEYSAM\n" +
//                 "HLA-A*11:01\tA59T+Q61R\tLLDILDTTGR\n" +
//                 "HLA-A*02:01\tQ61R+E62G\tILDTAGRGEY\n" +
//                 "HLA-A*32:01\tQ61R+E62G\tILDTAGRGEY\n" +
//                 "HLA-B*15:01\tA59T+Q61L\tTTGLEEYSAM\n" +
//                 "HLA-A*32:01\tA59T+E62G\tILDTTGQGEY\n" +
//                 "HLA-B*44:02\tQ61R+E62G\tILDTAGRGEY\n" +
//                 "HLA-A*23:01\tA59T+Q61H\tILDTTGHEEY\n" +
//                 "HLA-B*53:01\tA59T+Q61R\tTGREEYSAM\n" +
//                 "HLA-B*44:02\tA59T+E62G\tILDTTGQGEY\n" +
//                 "HLA-B*53:01\tQ61R+E62G\tILDTAGRGEY\n" +
//                 "HLA-B*44:03\tA59T+E62G\tILDTTGQGEY\n" +
//                 "HLA-A*11:01\tQ61R+E62G\tILDTAGRGEY\n" +
//                 "HLA-B*07:02\tA59T+Q61H\tTTGHEEYSAM\n" +
//                 "HLA-B*35:01\tQ61R+E62G\tTAGRGEYSA\n" +
//                 "HLA-B*07:02\tA59T+Q61L\tLDILDTTGL\n" +
//                 "HLA-B*07:02\tA59T+Q61H\tILDTTGHEEY\n" +
//                 "HLA-A*26:01\tQ61L+E62G\tAGLGEYSAM\n" +
//                 "HLA-A*26:01\tA59T+Q61L\tLDILDTTGL\n" +
//                 "HLA-B*44:03\tQ61R+E62G\tILDTAGRGEY\n" +
//                 "HLA-B*44:02\tA59T+Q61L\tLDILDTTGL\n" +
//                 "HLA-A*24:02\tA59T+Q61H\tILDTTGHEEY\n" +
//                 "HLA-B*40:01\tA59T+Q61H\tTGHEEYSAM\n" +
//                 "HLA-B*53:01\tQ61H+E62G\tLDTAGHGEY\n" +
//                 "HLA-B*44:03\tA59T+Q61L\tLDILDTTGL\n" +
//                 "HLA-B*53:01\tA59T+Q61H\tTTGHEEYSAM\n" +
//                 "HLA-A*23:01\tA59T+Q61L\tTGLEEYSAM\n" +
//                 "HLA-A*24:02\tA59T+Q61L\tTGLEEYSAM\n" +
//                 "HLA-A*23:01\tA59T+E62G\tILDTTGQGEY\n" +
//                 "HLA-A*24:02\tA59T+Q61H\tTGHEEYSAM\n";
        String data = args[0];

        String[] rowArray = data.split("\\r?\\n");
        String[] alleles = new String[rowArray.length];
        String[] epitopes = new String[rowArray.length];

        for(int i = 0; i < rowArray.length; i++)
        {
            alleles[i] = rowArray[i].substring(0,11);
            epitopes[i] = rowArray[i].substring(rowArray[i].lastIndexOf("\t")+1);
        }

        String[] epitopesNR = new HashSet<String>(Arrays.asList(epitopes)).toArray(new String[0]);
        String[][] alleleList = new String[epitopesNR.length][0];

        for(int i = 0; i < epitopesNR.length; i++)
        {
            String tempEpi = epitopesNR[i];

            for(int w = 0; w < epitopes.length; w++)
            {
                if(epitopes[w].equals(tempEpi))
                {
                    String[] tempArray = new String[alleleList[i].length+1];
                    for(int k = 0; k < alleleList[i].length; k++)
                    {
                        tempArray[k] = alleleList[i][k];
                    }
                    tempArray[tempArray.length-1] = alleles[w];
                    alleleList[i] = tempArray;
                }
            }
        }

        ArrayList<String> optimEpi = new ArrayList<String>();
        HashSet<String> allelesNRARRAY = new HashSet<String>(Arrays.asList(alleles));
        ArrayList<String> allelesNR = new ArrayList<String>();
        for(String s:allelesNRARRAY)
        {
            if(!allelesNR.contains(s))
            {
                allelesNR.add(s);
            }
        }
        ArrayList <String[]> optimAllele = new ArrayList<String[]>();

        String searchAllele;
        ArrayList<Integer> matchIndex = new ArrayList<Integer>();
        int index = 0;
        int temp = 0;

        while(allelesNR.size() > 0)
        {
            searchAllele = allelesNR.get(0);
            for(int i = 0; i < alleleList.length; i++)
            {
                for(int w = 0; w < alleleList[i].length; w++)
                {
                    if(alleleList[i][w].equals(searchAllele))
                    {
                        matchIndex.add(i);
                    }
                }
            }
            for(int i = 0; i < matchIndex.size(); i++)
            {
                if(alleleList[matchIndex.get(i)].length > temp)
                {
                    temp = alleleList[matchIndex.get(i)].length;
                    index = matchIndex.get(i);
                }
            }
            optimEpi.add(epitopesNR[index]);
            optimAllele.add(alleleList[index]);

			/*for(String[] s1:optimAllele)
			{
				for(String s2:s1)
				{
					System.out.print(s2 + " ");
				}
				System.out.println();
			}*/

            for(String s:alleleList[index])
            {
                allelesNR.remove(s);
            }
            matchIndex.clear();
            index = 0;
            temp = 0;
        }

        for(int i = 0; i < optimEpi.size(); i++)
        {
            System.out.print(optimEpi.get(i) + " ");
            for(int w = 0; w < optimAllele.get(i).length-1; w++)
            {
                System.out.print(optimAllele.get(i)[w] + ",");
            }
            System.out.print(optimAllele.get(i)[optimAllele.get(i).length-1]);
            System.out.println();
        }
    }

}