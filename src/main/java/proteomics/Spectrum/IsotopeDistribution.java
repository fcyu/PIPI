package proteomics.Spectrum;

import java.util.*;

class IsotopeDistribution {

    private static final double averagineC = 4.9384;
    private static final double averagineH = 7.7583;
    private static final double averagineN = 1.3577;
    private static final double averagineO = 1.4773;
    private static final double averagineS = 0.0417;

    private Map<String, Peak[]> elementIsotopeMap = new HashMap<>();
    private final double limit;
    private final double inverseAveragineMonoMass;
    private final Map<String, Double> elementTable;

    IsotopeDistribution(Map<String, Double> elementTable, double limit, String labelling) {
        this.elementTable = elementTable;
        this.limit = limit;
        inverseAveragineMonoMass = 1 / (averagineC * elementTable.get("C") + averagineH * elementTable.get("H") + averagineN * elementTable.get("N") + averagineO * elementTable.get("O") + averagineS * elementTable.get("S"));

        Peak[] peakArray = new Peak[2];
        peakArray[0] = new Peak(1.0078246, 0.99985);
        peakArray[1] = new Peak(2.0141021, 0.00015);
        elementIsotopeMap.put("H", peakArray);

        peakArray = new Peak[2];
        peakArray[0] = new Peak(3.01603, 0.00000138);
        peakArray[1] = new Peak(4.00260, 0.99999862);
        elementIsotopeMap.put("He", peakArray);

        peakArray = new Peak[2];
        peakArray[0] = new Peak(6.015121, 0.075);
        peakArray[1] = new Peak(7.016003, 0.925);
        elementIsotopeMap.put("Li", peakArray);

        peakArray = new Peak[1];
        peakArray[0] = new Peak(9.012182, 1.0);
        elementIsotopeMap.put("Be", peakArray);

        peakArray = new Peak[2];
        peakArray[0] = new Peak(10.012937, 0.199);
        peakArray[1] = new Peak(11.009305, 0.801);
        elementIsotopeMap.put("B", peakArray);

        peakArray = new Peak[2];
        peakArray[0] = new Peak(12.0000000, 0.988930);
        peakArray[1] = new Peak(13.0033554, 0.011070);
        elementIsotopeMap.put("C", peakArray);

        if (labelling.contentEquals("N15")) {
            peakArray = new Peak[1];
            peakArray[0] = new Peak(15.0001088, 1);
            elementIsotopeMap.put("N", peakArray);
        } else {
            peakArray = new Peak[2];
            peakArray[0] = new Peak(14.0030732, 0.996337);
            peakArray[1] = new Peak(15.0001088, 0.003663);
            elementIsotopeMap.put("N", peakArray);
        }

        peakArray = new Peak[3];
        peakArray[0] = new Peak(15.9949141, 0.997590);
        peakArray[1] = new Peak(16.9991322, 0.000374);
        peakArray[2] = new Peak(17.9991616, 0.002036);
        elementIsotopeMap.put("O", peakArray);

        peakArray = new Peak[1];
        peakArray[0] = new Peak(18.9984032, 1.0);
        elementIsotopeMap.put("F", peakArray);

        peakArray = new Peak[3];
        peakArray[0] = new Peak(19.992435, 0.9048);
        peakArray[1] = new Peak(20.993843, 0.0027);
        peakArray[2] = new Peak(21.991383, 0.0925);
        elementIsotopeMap.put("Ne", peakArray);

        peakArray = new Peak[1];
        peakArray[0] = new Peak(22.989767, 1.0);
        elementIsotopeMap.put("Na", peakArray);

        peakArray = new Peak[3];
        peakArray[0] = new Peak(23.985042, 0.7899);
        peakArray[1] = new Peak(24.985837, 0.1000);
        peakArray[2] = new Peak(25.982593, 0.1101);
        elementIsotopeMap.put("Mg", peakArray);

        peakArray = new Peak[1];
        peakArray[0] = new Peak(26.981539, 1.0);
        elementIsotopeMap.put("Al", peakArray);

        peakArray = new Peak[3];
        peakArray[0] = new Peak(27.976927, 0.9223);
        peakArray[1] = new Peak(28.976495, 0.0467);
        peakArray[2] = new Peak(29.973770, 0.0310);
        elementIsotopeMap.put("Si", peakArray);

        peakArray = new Peak[1];
        peakArray[0] = new Peak(30.973762, 1.0);
        elementIsotopeMap.put("P", peakArray);

        peakArray = new Peak[4];
        peakArray[0] = new Peak(31.972070, 0.9502);
        peakArray[1] = new Peak(32.971456, 0.0075);
        peakArray[2] = new Peak(33.967866, 0.0421);
        peakArray[3] = new Peak(35.967080, 0.0002);
        elementIsotopeMap.put("S", peakArray);

        peakArray = new Peak[2];
        peakArray[0] = new Peak(34.9688531, 0.755290);
        peakArray[1] = new Peak(36.9659034, 0.244710);
        elementIsotopeMap.put("Cl", peakArray);

        peakArray = new Peak[3];
        peakArray[0] = new Peak(35.967545, 0.00337);
        peakArray[1] = new Peak(37.962732, 0.00063);
        peakArray[2] = new Peak(39.962384, 0.99600);
        elementIsotopeMap.put("Ar", peakArray);

        peakArray = new Peak[3];
        peakArray[0] = new Peak(38.963707, 0.932581);
        peakArray[1] = new Peak(39.963999, 0.000117);
        peakArray[2] = new Peak(40.961825, 0.067302);
        elementIsotopeMap.put("K", peakArray);

        peakArray = new Peak[6];
        peakArray[0] = new Peak(39.962591, 0.96941);
        peakArray[1] = new Peak(41.958618, 0.00647);
        peakArray[2] = new Peak(42.958766, 0.00135);
        peakArray[3] = new Peak(43.955480, 0.02086);
        peakArray[4] = new Peak(45.953689, 0.00004);
        peakArray[5] = new Peak(47.952533, 0.00187);
        elementIsotopeMap.put("Ca", peakArray);

        peakArray = new Peak[1];
        peakArray[0] = new Peak(44.955910, 1.0);
        elementIsotopeMap.put("Sc", peakArray);

        peakArray = new Peak[5];
        peakArray[0] = new Peak(45.952629, 0.080);
        peakArray[1] = new Peak(46.951764, 0.073);
        peakArray[2] = new Peak(47.947947, 0.738);
        peakArray[3] = new Peak(48.947871, 0.055);
        peakArray[4] = new Peak(49.944792, 0.054);
        elementIsotopeMap.put("Ti", peakArray);

        peakArray = new Peak[2];
        peakArray[0] = new Peak(49.947161, 0.00250);
        peakArray[1] = new Peak(50.943962, 0.99750);
        elementIsotopeMap.put("V", peakArray);

        peakArray = new Peak[4];
        peakArray[0] = new Peak(49.946046, 0.04345);
        peakArray[1] = new Peak(51.940509, 0.83790);
        peakArray[2] = new Peak(52.940651, 0.09500);
        peakArray[3] = new Peak(53.938882, 0.02365);
        elementIsotopeMap.put("Cr", peakArray);

        peakArray = new Peak[1];
        peakArray[0] = new Peak(54.938047, 1.0);
        elementIsotopeMap.put("Mn", peakArray);

        peakArray = new Peak[4];
        peakArray[0] = new Peak(53.939612, 0.0590);
        peakArray[1] = new Peak(55.934939, 0.9172);
        peakArray[2] = new Peak(56.935396, 0.0210);
        peakArray[3] = new Peak(57.933277, 0.0028);
        elementIsotopeMap.put("Fe", peakArray);

        peakArray = new Peak[1];
        peakArray[0] = new Peak(58.933198, 1.0);
        elementIsotopeMap.put("Co", peakArray);

        peakArray = new Peak[5];
        peakArray[0] = new Peak(57.935346, 0.6827);
        peakArray[1] = new Peak(59.930788, 0.2610);
        peakArray[2] = new Peak(60.931058, 0.0113);
        peakArray[3] = new Peak(61.928346, 0.0359);
        peakArray[4] = new Peak(63.927968, 0.0091);
        elementIsotopeMap.put("Ni", peakArray);

        peakArray = new Peak[2];
        peakArray[0] = new Peak(62.939598, 0.6917);
        peakArray[1] = new Peak(64.927793, 0.3083);
        elementIsotopeMap.put("Cu", peakArray);

        peakArray = new Peak[5];
        peakArray[0] = new Peak(63.929145, 0.486);
        peakArray[1] = new Peak(65.926034, 0.279);
        peakArray[2] = new Peak(66.927129, 0.041);
        peakArray[3] = new Peak(67.924846, 0.188);
        peakArray[4] = new Peak(69.925325, 0.006);
        elementIsotopeMap.put("Zn", peakArray);

        peakArray = new Peak[2];
        peakArray[0] = new Peak(68.925580, 0.60108);
        peakArray[1] = new Peak(70.924700, 0.39892);
        elementIsotopeMap.put("Ga", peakArray);

        peakArray = new Peak[5];
        peakArray[0] = new Peak(69.924250, 0.205);
        peakArray[1] = new Peak(71.922079, 0.274);
        peakArray[2] = new Peak(72.923463, 0.078);
        peakArray[3] = new Peak(73.921177, 0.365);
        peakArray[4] = new Peak(75.921401, 0.078);
        elementIsotopeMap.put("Ge", peakArray);

        peakArray = new Peak[1];
        peakArray[0] = new Peak(74.921594, 1.0);
        elementIsotopeMap.put("As", peakArray);

        peakArray = new Peak[6];
        peakArray[0] = new Peak(73.922475, 0.009);
        peakArray[1] = new Peak(75.919212, 0.091);
        peakArray[2] = new Peak(76.919912, 0.076);
        peakArray[3] = new Peak(77.9190, 0.236);
        peakArray[4] = new Peak(79.916520, 0.499);
        peakArray[5] = new Peak(81.916698, 0.089);
        elementIsotopeMap.put("Se", peakArray);

        peakArray = new Peak[2];
        peakArray[0] = new Peak(78.918336, 0.5069);
        peakArray[1] = new Peak(80.916289, 0.4931);
        elementIsotopeMap.put("Br", peakArray);

        peakArray = new Peak[6];
        peakArray[0] = new Peak(77.914, 0.0035);
        peakArray[1] = new Peak(79.916380, 0.0225);
        peakArray[2] = new Peak(81.913482, 0.116);
        peakArray[3] = new Peak(82.914135, 0.115);
        peakArray[4] = new Peak(83.911507, 0.570);
        peakArray[5] = new Peak(85.910616, 0.173);
        elementIsotopeMap.put("Kr", peakArray);

        peakArray = new Peak[2];
        peakArray[0] = new Peak(84.911794, 0.7217);
        peakArray[1] = new Peak(86.909187, 0.2783);
        elementIsotopeMap.put("Rb", peakArray);

        peakArray = new Peak[4];
        peakArray[0] = new Peak(83.913430, 0.0056);
        peakArray[1] = new Peak(85.909267, 0.0986);
        peakArray[2] = new Peak(86.908884, 0.0700);
        peakArray[3] = new Peak(87.905619, 0.8258);
        elementIsotopeMap.put("Sr", peakArray);

        peakArray = new Peak[1];
        peakArray[0] = new Peak(88.905849, 1.0);
        elementIsotopeMap.put("Y", peakArray);

        peakArray = new Peak[5];
        peakArray[0] = new Peak(89.904703, 0.5145);
        peakArray[1] = new Peak(90.905644, 0.1122);
        peakArray[2] = new Peak(91.905039, 0.1715);
        peakArray[3] = new Peak(93.906314, 0.1738);
        peakArray[4] = new Peak(95.908275, 0.0280);
        elementIsotopeMap.put("Zr", peakArray);

        peakArray = new Peak[1];
        peakArray[0] = new Peak(92.906377, 1.0);
        elementIsotopeMap.put("Nb", peakArray);

        peakArray = new Peak[7];
        peakArray[0] = new Peak(91.906808, 0.1484);
        peakArray[1] = new Peak(93.905085, 0.0925);
        peakArray[2] = new Peak(94.905840, 0.1592);
        peakArray[3] = new Peak(95.904678, 0.1668);
        peakArray[4] = new Peak(96.906020, 0.0955);
        peakArray[5] = new Peak(97.905406, 0.2413);
        peakArray[6] = new Peak(99.907477, 0.0963);
        elementIsotopeMap.put("Mo", peakArray);

        peakArray = new Peak[1];
        peakArray[0] = new Peak(98.0, 1.0);
        elementIsotopeMap.put("Tc", peakArray);

        peakArray = new Peak[7];
        peakArray[0] = new Peak(95.907599, 0.0554);
        peakArray[1] = new Peak(97.905287, 0.0186);
        peakArray[2] = new Peak(98.905939, 0.127);
        peakArray[3] = new Peak(99.904219, 0.126);
        peakArray[4] = new Peak(100.905582, 0.171);
        peakArray[5] = new Peak(101.904348, 0.316);
        peakArray[6] = new Peak(103.905424, 0.186);
        elementIsotopeMap.put("Ru", peakArray);

        peakArray = new Peak[1];
        peakArray[0] = new Peak(102.905500, 1.0);
        elementIsotopeMap.put("Rh", peakArray);

        peakArray = new Peak[6];
        peakArray[0] = new Peak(101.905634, 0.0102);
        peakArray[1] = new Peak(103.904029, 0.1114);
        peakArray[2] = new Peak(104.905079, 0.2233);
        peakArray[3] = new Peak(105.903478, 0.2733);
        peakArray[4] = new Peak(107.903895, 0.2646);
        peakArray[5] = new Peak(109.905167, 0.1172);
        elementIsotopeMap.put("Pd", peakArray);

        peakArray = new Peak[2];
        peakArray[0] = new Peak(106.905092, 0.51839);
        peakArray[1] = new Peak(108.904757, 0.48161);
        elementIsotopeMap.put("Ag", peakArray);

        peakArray = new Peak[8];
        peakArray[0] = new Peak(105.906461, 0.0125);
        peakArray[1] = new Peak(107.904176, 0.0089);
        peakArray[2] = new Peak(109.903005, 0.1249);
        peakArray[3] = new Peak(110.904182, 0.1280);
        peakArray[4] = new Peak(111.902758, 0.2413);
        peakArray[5] = new Peak(112.904400, 0.1222);
        peakArray[6] = new Peak(113.903357, 0.2873);
        peakArray[7] = new Peak(115.904754, 0.0749);
        elementIsotopeMap.put("Cd", peakArray);

        peakArray = new Peak[2];
        peakArray[0] = new Peak(112.904061, 0.043);
        peakArray[1] = new Peak(114.903880, 0.957);
        elementIsotopeMap.put("In", peakArray);

        peakArray = new Peak[10];
        peakArray[0] = new Peak(111.904826, 0.0097);
        peakArray[1] = new Peak(113.902784, 0.0065);
        peakArray[2] = new Peak(114.903348, 0.0036);
        peakArray[3] = new Peak(115.901747, 0.1453);
        peakArray[4] = new Peak(116.902956, 0.0768);
        peakArray[5] = new Peak(117.901609, 0.2422);
        peakArray[6] = new Peak(118.903310, 0.0858);
        peakArray[7] = new Peak(119.902200, 0.3259);
        peakArray[8] = new Peak(121.903440, 0.0463);
        peakArray[9] = new Peak(123.905274, 0.0579);
        elementIsotopeMap.put("Sn", peakArray);

        peakArray = new Peak[2];
        peakArray[0] = new Peak(120.903821, 0.574);
        peakArray[1] = new Peak(122.904216, 0.426);
        elementIsotopeMap.put("Sb", peakArray);

        peakArray = new Peak[8];
        peakArray[0] = new Peak(119.904048, 0.00095);
        peakArray[1] = new Peak(121.903054, 0.0259);
        peakArray[2] = new Peak(122.904271, 0.00905);
        peakArray[3] = new Peak(123.902823, 0.0479);
        peakArray[4] = new Peak(124.904433, 0.0712);
        peakArray[5] = new Peak(125.903314, 0.1893);
        peakArray[6] = new Peak(127.904463, 0.3170);
        peakArray[7] = new Peak(129.906229, 0.3387);
        elementIsotopeMap.put("Te", peakArray);

        peakArray = new Peak[1];
        peakArray[0] = new Peak(126.904473, 1.0);
        elementIsotopeMap.put("I", peakArray);

        peakArray = new Peak[9];
        peakArray[0] = new Peak(123.905894, 0.0010);
        peakArray[1] = new Peak(125.904281, 0.0009);
        peakArray[2] = new Peak(127.903531, 0.0191);
        peakArray[3] = new Peak(128.904780, 0.264);
        peakArray[4] = new Peak(129.903509, 0.041);
        peakArray[5] = new Peak(130.905072, 0.212);
        peakArray[6] = new Peak(131.904144, 0.269);
        peakArray[7] = new Peak(133.905395, 0.104);
        peakArray[8] = new Peak(135.907214, 0.089);
        elementIsotopeMap.put("Xe", peakArray);

        peakArray = new Peak[1];
        peakArray[0] = new Peak(132.905429, 1.0);
        elementIsotopeMap.put("Cs", peakArray);

        peakArray = new Peak[7];
        peakArray[0] = new Peak(129.906282, 0.00106);
        peakArray[1] = new Peak(131.905042, 0.00101);
        peakArray[2] = new Peak(133.904486, 0.0242);
        peakArray[3] = new Peak(134.905665, 0.06593);
        peakArray[4] = new Peak(135.904553, 0.0785);
        peakArray[5] = new Peak(136.905812, 0.1123);
        peakArray[6] = new Peak(137.905232, 0.7170);
        elementIsotopeMap.put("Ba", peakArray);

        peakArray = new Peak[2];
        peakArray[0] = new Peak(137.90711, 0.00090);
        peakArray[1] = new Peak(138.906347, 0.99910);
        elementIsotopeMap.put("La", peakArray);

        peakArray = new Peak[4];
        peakArray[0] = new Peak(135.907140, 0.0019);
        peakArray[1] = new Peak(137.905985, 0.0025);
        peakArray[2] = new Peak(139.905433, 0.8843);
        peakArray[3] = new Peak(141.909241, 0.1113);
        elementIsotopeMap.put("Ce", peakArray);

        peakArray = new Peak[1];
        peakArray[0] = new Peak(140.907647, 1.0);
        elementIsotopeMap.put("Pr", peakArray);

        peakArray = new Peak[7];
        peakArray[0] = new Peak(141.907719, 0.2713);
        peakArray[1] = new Peak(142.909810, 0.1218);
        peakArray[2] = new Peak(143.910083, 0.2380);
        peakArray[3] = new Peak(144.912570, 0.0830);
        peakArray[4] = new Peak(145.913113, 0.1719);
        peakArray[5] = new Peak(147.916889, 0.0576);
        peakArray[6] = new Peak(149.920887, 0.0564);
        elementIsotopeMap.put("Nd", peakArray);

        peakArray = new Peak[1];
        peakArray[0] = new Peak(145.0, 1.0);
        elementIsotopeMap.put("Pm", peakArray);

        peakArray = new Peak[7];
        peakArray[0] = new Peak(143.911998, 0.031);
        peakArray[1] = new Peak(146.914895, 0.150);
        peakArray[2] = new Peak(147.914820, 0.113);
        peakArray[3] = new Peak(148.917181, 0.138);
        peakArray[4] = new Peak(149.917273, 0.074);
        peakArray[5] = new Peak(151.919729, 0.267);
        peakArray[6] = new Peak(153.922206, 0.227);
        elementIsotopeMap.put("Sm", peakArray);

        peakArray = new Peak[2];
        peakArray[0] = new Peak(150.919847, 0.478);
        peakArray[1] = new Peak(152.921225, 0.522);
        elementIsotopeMap.put("Eu", peakArray);

        peakArray = new Peak[7];
        peakArray[0] = new Peak(151.919786, 0.0020);
        peakArray[1] = new Peak(153.920861, 0.0218);
        peakArray[2] = new Peak(154.922618, 0.1480);
        peakArray[3] = new Peak(155.922118, 0.2047);
        peakArray[4] = new Peak(156.923956, 0.1565);
        peakArray[5] = new Peak(157.924099, 0.2484);
        peakArray[6] = new Peak(159.927049, 0.2186);
        elementIsotopeMap.put("Gd", peakArray);

        peakArray = new Peak[1];
        peakArray[0] = new Peak(158.925342, 1.0);
        elementIsotopeMap.put("Tb", peakArray);

        peakArray = new Peak[7];
        peakArray[0] = new Peak(155.925277, 0.0006);
        peakArray[1] = new Peak(157.924403, 0.0010);
        peakArray[2] = new Peak(159.925193, 0.0234);
        peakArray[3] = new Peak(160.926930, 0.189);
        peakArray[4] = new Peak(161.926795, 0.255);
        peakArray[5] = new Peak(162.928728, 0.249);
        peakArray[6] = new Peak(163.929171, 0.282);
        elementIsotopeMap.put("Dy", peakArray);

        peakArray = new Peak[1];
        peakArray[0] = new Peak(164.930319, 1.0);
        elementIsotopeMap.put("Ho", peakArray);

        peakArray = new Peak[6];
        peakArray[0] = new Peak(161.928775, 0.0014);
        peakArray[1] = new Peak(163.929198, 0.0161);
        peakArray[2] = new Peak(165.930290, 0.336);
        peakArray[3] = new Peak(166.932046, 0.2295);
        peakArray[4] = new Peak(167.932368, 0.268);
        peakArray[5] = new Peak(169.935461, 0.149);
        elementIsotopeMap.put("Er", peakArray);

        peakArray = new Peak[1];
        peakArray[0] = new Peak(168.934212, 1.0);
        elementIsotopeMap.put("Tm", peakArray);

        peakArray = new Peak[7];
        peakArray[0] = new Peak(167.933894, 0.0013);
        peakArray[1] = new Peak(169.934759, 0.0305);
        peakArray[2] = new Peak(170.936323, 0.143);
        peakArray[3] = new Peak(171.936378, 0.219);
        peakArray[4] = new Peak(172.938208, 0.1612);
        peakArray[5] = new Peak(173.938859, 0.318);
        peakArray[6] = new Peak(175.942564, 0.127);
        elementIsotopeMap.put("Yb", peakArray);

        peakArray = new Peak[2];
        peakArray[0] = new Peak(174.940770, 0.9741);
        peakArray[1] = new Peak(175.942679, 0.0259);
        elementIsotopeMap.put("Lu", peakArray);

        peakArray = new Peak[6];
        peakArray[0] = new Peak(173.940044, 0.00162);
        peakArray[1] = new Peak(175.941406, 0.05206);
        peakArray[2] = new Peak(176.943217, 0.18606);
        peakArray[3] = new Peak(177.943696, 0.27297);
        peakArray[4] = new Peak(178.945812, 0.13629);
        peakArray[5] = new Peak(179.946545, 0.35100);
        elementIsotopeMap.put("Hf", peakArray);

        peakArray = new Peak[2];
        peakArray[0] = new Peak(179.947462, 0.00012);
        peakArray[1] = new Peak(180.947992, 0.99988);
        elementIsotopeMap.put("Ta", peakArray);

        peakArray = new Peak[5];
        peakArray[0] = new Peak(179.946701, 0.0012);
        peakArray[1] = new Peak(181.948202, 0.263);
        peakArray[2] = new Peak(182.950220, 0.1428);
        peakArray[3] = new Peak(183.950928, 0.307);
        peakArray[4] = new Peak(185.954357, 0.286);
        elementIsotopeMap.put("W", peakArray);

        peakArray = new Peak[2];
        peakArray[0] = new Peak(184.952951, 0.3740);
        peakArray[1] = new Peak(186.955744, 0.6260);
        elementIsotopeMap.put("Re", peakArray);

        peakArray = new Peak[7];
        peakArray[0] = new Peak(183.952488, 0.0002);
        peakArray[1] = new Peak(185.953830, 0.0158);
        peakArray[2] = new Peak(186.955741, 0.016);
        peakArray[3] = new Peak(187.955860, 0.133);
        peakArray[4] = new Peak(188.958137, 0.161);
        peakArray[5] = new Peak(189.958436, 0.264);
        peakArray[6] = new Peak(191.961467, 0.410);
        elementIsotopeMap.put("Os", peakArray);

        peakArray = new Peak[2];
        peakArray[0] = new Peak(190.960584, 0.373);
        peakArray[1] = new Peak(192.962917, 0.627);
        elementIsotopeMap.put("Ir", peakArray);

        peakArray = new Peak[6];
        peakArray[0] = new Peak(189.959917, 0.0001);
        peakArray[1] = new Peak(191.961019, 0.0079);
        peakArray[2] = new Peak(193.962655, 0.329);
        peakArray[3] = new Peak(194.964766, 0.338);
        peakArray[4] = new Peak(195.964926, 0.253);
        peakArray[5] = new Peak(197.967869, 0.072);
        elementIsotopeMap.put("Pt", peakArray);

        peakArray = new Peak[1];
        peakArray[0] = new Peak(196.966543, 1.0);
        elementIsotopeMap.put("Au", peakArray);

        peakArray = new Peak[7];
        peakArray[0] = new Peak(195.965807, 0.0015);
        peakArray[1] = new Peak(197.966743, 0.100);
        peakArray[2] = new Peak(198.968254, 0.169);
        peakArray[3] = new Peak(199.968300, 0.231);
        peakArray[4] = new Peak(200.970277, 0.132);
        peakArray[5] = new Peak(201.970617, 0.298);
        peakArray[6] = new Peak(203.973467, 0.0685);
        elementIsotopeMap.put("Hg", peakArray);

        peakArray = new Peak[2];
        peakArray[0] = new Peak(202.972320, 0.29524);
        peakArray[1] = new Peak(204.974401, 0.70476);
        elementIsotopeMap.put("Tl", peakArray);

        peakArray = new Peak[4];
        peakArray[0] = new Peak(203.973020, 0.014);
        peakArray[1] = new Peak(205.974440, 0.241);
        peakArray[2] = new Peak(206.975872, 0.221);
        peakArray[3] = new Peak(207.976627, 0.524);
        elementIsotopeMap.put("Pb", peakArray);

        peakArray = new Peak[1];
        peakArray[0] = new Peak(208.980374, 1.0);
        elementIsotopeMap.put("Bi", peakArray);

        peakArray = new Peak[1];
        peakArray[0] = new Peak(209.0, 1.0);
        elementIsotopeMap.put("Po", peakArray);

        peakArray = new Peak[1];
        peakArray[0] = new Peak(210.0, 1.0);
        elementIsotopeMap.put("At", peakArray);

        peakArray = new Peak[1];
        peakArray[0] = new Peak(222.0, 1.0);
        elementIsotopeMap.put("Rn", peakArray);

        peakArray = new Peak[1];
        peakArray[0] = new Peak(223.0, 1.0);
        elementIsotopeMap.put("Fr", peakArray);

        peakArray = new Peak[1];
        peakArray[0] = new Peak(226.025, 1.0);
        elementIsotopeMap.put("Ra", peakArray);

        peakArray = new Peak[1];
        peakArray[0] = new Peak(227.028, 1.0);
        elementIsotopeMap.put("Ac", peakArray);

        peakArray = new Peak[1];
        peakArray[0] = new Peak(232.038054, 1.0);
        elementIsotopeMap.put("Th", peakArray);

        peakArray = new Peak[1];
        peakArray[0] = new Peak(231.0359, 1.0);
        elementIsotopeMap.put("Pa", peakArray);

        peakArray = new Peak[3];
        peakArray[0] = new Peak(234.040946, 0.000055);
        peakArray[1] = new Peak(235.043924, 0.00720);
        peakArray[2] = new Peak(238.050784, 0.992745);
        elementIsotopeMap.put("U", peakArray);

        peakArray = new Peak[1];
        peakArray[0] = new Peak(237.048, 1.0);
        elementIsotopeMap.put("Np", peakArray);

        peakArray = new Peak[1];
        peakArray[0] = new Peak(244.0, 1.0);
        elementIsotopeMap.put("Pu", peakArray);

        peakArray = new Peak[1];
        peakArray[0] = new Peak(243.0, 1.0);
        elementIsotopeMap.put("Am", peakArray);

        peakArray = new Peak[1];
        peakArray[0] = new Peak(247.0, 1.0);
        elementIsotopeMap.put("Cm", peakArray);

        peakArray = new Peak[1];
        peakArray[0] = new Peak(247.0, 1.0);
        elementIsotopeMap.put("Bk", peakArray);

        peakArray = new Peak[1];
        peakArray[0] = new Peak(251.0, 1.0);
        elementIsotopeMap.put("Cf", peakArray);

        peakArray = new Peak[1];
        peakArray[0] = new Peak(252.0, 1.0);
        elementIsotopeMap.put("Es", peakArray);

        peakArray = new Peak[1];
        peakArray[0] = new Peak(257.0, 1.0);
        elementIsotopeMap.put("Fm", peakArray);

        peakArray = new Peak[1];
        peakArray[0] = new Peak(258.0, 1.0);
        elementIsotopeMap.put("Md", peakArray);

        peakArray = new Peak[1];
        peakArray[0] = new Peak(259.0, 1.0);
        elementIsotopeMap.put("No", peakArray);

        peakArray = new Peak[1];
        peakArray[0] = new Peak(260.0, 1.0);
        elementIsotopeMap.put("Lr", peakArray);
    }

    Map<String, Integer> getElementMapFromMonoMass(double mass) {
        double unit = mass * inverseAveragineMonoMass;
        int C = (int) Math.round(averagineC * unit);
        int N = (int) Math.round(averagineN * unit);
        int O = (int) Math.round(averagineO * unit);
        int S = (int) Math.round(averagineS * unit);
        int H = (int) Math.round((mass - C * elementTable.get("C") - N * elementTable.get("N") - O * elementTable.get("O") - S * elementTable.get("S")) / elementTable.get("H"));

        Map<String, Integer> elementMap = new HashMap<>();
        if (C > 0) {
            elementMap.put("C", C);
        }
        if (N > 0) {
            elementMap.put("N", N);
        }
        if (O > 0) {
            elementMap.put("O", O);
        }
        if (S > 0) {
            elementMap.put("S", S);
        }
        if (H > 0) {
            elementMap.put("H", H);
        }
        return elementMap;
    }

    List<Peak> calculate(Map<String, Integer> formMap) {
        List<Peak> result = new ArrayList<>();
        result.add(new Peak(0, 1));
        for (String element : formMap.keySet()) {
            List<List<Peak>> sal = new ArrayList<>();
            sal.add(Arrays.asList(elementIsotopeMap.get(element)));
            int n = formMap.get(element);
            int j = 0;
            while (n > 0) {
                if (j == sal.size()) {
                    sal.add(convolute(sal.get(j - 1), sal.get(j - 1)));
                    prune(sal.get(j), limit);
                }
                if ((n & 1) != 0) {
                    result = convolute(result, sal.get(j));
                    prune(result, limit);
                }
                n >>= 1;
                ++j;
            }
        }

        prune(result, 1e-6); // prune the final result to a reasonable precision.
        return result;
    }

    private List<Peak> convolute(List<Peak> g, List<Peak> f) {
        List<Peak> h = new ArrayList<>();
        int gN = g.size();
        int fN = f.size();
        if (gN != 0 && fN != 0) {
            for (int k = 0; k < gN + fN - 1; ++k) {
                double sumWeight = 0;
                double sumMass = 0;
                int start = Math.max(0, k - fN + 1);
                int end = Math.min(gN - 1, k);
                for (int i = start; i <= end; ++i) {
                    double weight = g.get(i).realArea * f.get(k - i).realArea;
                    double mass = g.get(i).mass + f.get(k - i).mass;
                    sumWeight += weight;
                    sumMass += weight * mass;
                }
                Peak peak;
                if (sumWeight == 0) {
                    peak = new Peak(-1, sumWeight);
                } else {
                    peak = new Peak(sumMass / sumWeight, sumWeight);
                }
                h.add(peak);
            }
        }
        return h;
    }

    private void prune(List<Peak> f, double limit) {
        List<Peak> toBeDeleteList = new LinkedList<>();
        for (Peak peak : f) {
            if (peak.realArea <= limit) {
                toBeDeleteList.add(peak);
            }
        }
        f.removeAll(toBeDeleteList);
    }


    public static class Peak {

        public final double mass;
        public final double realArea;

        public Peak(double mass, double realArea) {
            this.mass = mass;
            this.realArea = realArea;
        }
    }
}
