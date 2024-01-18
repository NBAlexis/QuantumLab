using System;
using System.Collections.Generic;
using System.IO;
using System.Runtime.Remoting.Contexts;
using System.Text.RegularExpressions;
using System.Xml;

namespace CMakeWrite
{
    class Program
    {
        enum EArch
        {
            //970
            EArcSM52,

            //RTX1080, 1070, 1060
            EArcSM61,

            //V100
            EArcSM70,

            //RTX2080, 2070, 2060, 1660
            EArcSM75,

            //RTX3060, 3070, 3080
            EArcSM86,
        }

        readonly static string[] ArchNames =
        {
            //GTX 970M
            "compute_52",
            "compute_61",
            "compute_70",
            "compute_75",
            "compute_86",
        };

        readonly static string[] CodeNames =
        {
            "sm_52",
            "sm_61",
            "sm_70",
            "sm_75",
            "sm_86",
        };

        class CProjFile
        {
            public CProjFile(string sName)
            {
                bool bProjFound = false;
                m_sName = sName;
                DirectoryInfo codeFolder = new DirectoryInfo(Path.Combine(new[] { System.AppDomain.CurrentDomain.BaseDirectory, "../../" }));
                FileInfo[] allProjects = codeFolder.GetFiles("*.vcxproj", SearchOption.AllDirectories);
                foreach (FileInfo projFiles in allProjects)
                {
                    Console.WriteLine(projFiles.Name);
                    if (projFiles.Name == sName + ".vcxproj")
                    {
                        bProjFound = true;
                        m_sProjectDir = projFiles.DirectoryName + "/";
                        m_sContent = File.ReadAllText(projFiles.FullName);
                    }
                }

                if (!bProjFound)
                {
                    Console.WriteLine(string.Format("Project {0} not found!\n", sName));
                    m_bExist = false;
                    return;
                }

                m_bExist = true;
                XmlDocument doc = new XmlDocument();
                doc.LoadXml(m_sContent);

                XmlNodeList allClInclude = doc.GetElementsByTagName("ClInclude");
                for (int i = 0; i < allClInclude.Count; ++i)
                {
                    m_lstAllHeaderFiles.Add(allClInclude[i].Attributes.GetNamedItem("Include").InnerText);
                }

                XmlNodeList allClCompile = doc.GetElementsByTagName("ClCompile");
                for (int i = 0; i < allClCompile.Count; ++i)
                {
                    if (null == allClCompile[i].Attributes.GetNamedItem("Include"))
                    {
                        //those are compile settings
                    }
                    else
                    {
                        m_lstAllCppFiles.Add(allClCompile[i].Attributes.GetNamedItem("Include").InnerText);
                    }
                }

                XmlNodeList allCuCompile = doc.GetElementsByTagName("CudaCompile");
                for (int i = 0; i < allCuCompile.Count; ++i)
                {
                    if (null == allCuCompile[i].Attributes.GetNamedItem("Include"))
                    {
                        //those are nvcc settings
                        ///<seealso cref="m_sContent">
                    }
                    else
                    {
                        m_lstAllCuFiles.Add(allCuCompile[i].Attributes.GetNamedItem("Include").InnerText);
                    }
                }
            }

            public int[] FindVersionNumber()
            {
                int[] ret = { 0, 0 };

                foreach (string sHeader in m_lstAllHeaderFiles)
                {
                    if (sHeader.Contains("CNDefine.h"))
                    {
                        string sFullPath = Path.GetFullPath(Path.Combine(m_sProjectDir, sHeader));
                        string sFileContent = File.ReadAllText(sFullPath);
                        Match allMatch = Regex.Match(sFileContent, @"__GVERSION[\s]+\(([\d])+\)");
                        if (allMatch.Success && allMatch.Groups.Count > 1)
                        {
                            int iMajor = 1;
                            int.TryParse(allMatch.Groups[1].ToString(), out iMajor);
                            ret[0] = iMajor;
                        }

                        allMatch = Regex.Match(sFileContent, @"__GVERSION_S[\s]+\(([\d])+\)");
                        if (allMatch.Success && allMatch.Groups.Count > 1)
                        {
                            int iMinor = 0;
                            int.TryParse(allMatch.Groups[1].ToString(), out iMinor);
                            ret[1] = iMinor;
                        }
                    }
                }

                return ret;
            }

            public readonly List<string> m_lstAllHeaderFiles = new List<string>();
            public readonly List<string> m_lstAllCuFiles = new List<string>();
            public readonly List<string> m_lstAllCppFiles = new List<string>();

            public string m_sName;
            public string m_sContent;
            public string m_sProjectDir;
            public bool m_bExist = false;

            public string GetCMakeContentAsApplication(string[] linkWith)
            {
                if (!m_bExist)
                {
                    return "";
                }
                string sAppName = m_sName;
                string sProjFolder = m_sProjectDir.Replace("\\", "/").Substring(0, m_sProjectDir.Length - 1);
                sProjFolder = Regex.Replace(sProjFolder, "([\\s\\S]*)Code/([\\s\\S]*)", "$2");

                var sRet = string.Format("\n\n\n# ==================== \n# {0} \n# =================\n\n", sAppName);
                if (m_sName.Equals("tests"))
                {
                    sRet += string.Format("include_directories({1}/{0}/tests)\n", sProjFolder, "${PROJECT_SOURCE_DIR}");
                }
                else
                {
                    sRet += string.Format("include_directories({1}/{0})\n", sProjFolder, "${PROJECT_SOURCE_DIR}");
                }
                    
                sRet += string.Format("add_executable({0} \n    ", sAppName);
                foreach (string sFileName in m_lstAllHeaderFiles)
                {
                    sRet += string.Format("{1}/{0}/" + sFileName.Replace("\\", "/") + "\n    ", sProjFolder, "${PROJECT_SOURCE_DIR}");
                }

                foreach (string sFileName in m_lstAllCppFiles)
                {
                    sRet += string.Format("{1}/{0}/" + sFileName.Replace("\\", "/") + "\n    ", sProjFolder, "${PROJECT_SOURCE_DIR}");
                }

                sRet += ")\n\n";

                sRet += string.Format("\ntarget_include_directories({0} PRIVATE {1})\n\n", m_sName, "${CMAKE_SOURCE_DIR}/../QuEST340/QuEST/include");
                if (m_sName.Equals("tests"))
                {
                    sRet += string.Format("\ntarget_include_directories({0} PRIVATE {1})\n\n", m_sName, "${CMAKE_SOURCE_DIR}/../QuEST340/tests/catch");
                }

                sRet += string.Format("target_compile_features({0} PUBLIC cxx_std_14)\n", sAppName);
                foreach (string linkWithName in linkWith)
                {
                    sRet += string.Format("target_link_libraries({0} {1})\n", sAppName, linkWithName);
                }
                return sRet;
            }

            public string GetCMakeContentAsCudaLib(string[] linkWith)
            {
                if (!m_bExist)
                {
                    return "";
                }
                string sProjFolder = m_sProjectDir.Replace("\\", "/").Substring(0, m_sProjectDir.Length - 1);
                sProjFolder = Regex.Replace(sProjFolder, "([\\s\\S]*)Code/([\\s\\S]*)", "$2");
                string sContent = string.Format("\n\n\n# ==================== \n# {0} \n# =================\n\n", m_sName);
                sContent += string.Format("include_directories(${1}/{0})\n", m_sName, "{PROJECT_SOURCE_DIR}");

                sContent += string.Format("add_library({0} STATIC\n    ", m_sName);
                foreach (string sFileName in m_lstAllHeaderFiles)
                {
                    sContent += string.Format("${2}/{0}/{1}\n    ", sProjFolder, sFileName.Replace("\\", "/"), "{PROJECT_SOURCE_DIR}");
                }
                foreach (string sFileName in m_lstAllCuFiles)
                {
                    sContent += string.Format("${2}/{0}/{1}\n    ", sProjFolder, sFileName.Replace("\\", "/"), "{PROJECT_SOURCE_DIR}");
                }
                foreach (string sFileName in m_lstAllCppFiles)
                {
                    sContent += string.Format("${2}/{0}/{1}\n    ", sProjFolder, sFileName.Replace("\\", "/"), "{PROJECT_SOURCE_DIR}");
                }
                sContent += ")\n\n";

                sContent += string.Format(@"# Request that {0} be built with -std=c++14
# As this is a public compile feature anything that links to 
# {0} will also build with -std=c++14
target_compile_features({0} PUBLIC cxx_std_14)
 
# We need to explicitly state that we need all CUDA files in the 
# {0} library to be built with -dc as the member functions 
# could be called by other libraries and executables
set_target_properties({0} PROPERTIES CUDA_SEPARABLE_COMPILATION ON)", m_sName);

                //Add include path here
                sContent += string.Format("\ntarget_include_directories({0} PRIVATE {1})\n\n", m_sName, "${CMAKE_SOURCE_DIR}/../QuEST340/QuEST/include");
                if (m_sName.Equals("QuEST"))
                {
                    sContent += string.Format("\ntarget_include_directories({0} PRIVATE {1})\n\n", m_sName, "${CMAKE_SOURCE_DIR}/../QuEST340/QuEST/src");
                }

                //Add links here
                sContent += string.Format("\n\ntarget_link_libraries({0} -lcurand)\n", m_sName);
                sContent += string.Format("target_link_libraries({0} -lcufft)\n", m_sName);
                sContent += string.Format("target_link_libraries({0} -lcusolver)\n", m_sName);
                sContent += string.Format("target_link_libraries({0} -lcublas)\n", m_sName);
                //sContent += string.Format("target_link_libraries({0} -lcudadevrt)\n", m_sName);
                foreach (string linkWithName in linkWith)
                {
                    sContent += string.Format("target_link_libraries({0} {1})\n", m_sName, linkWithName);
                }

                sContent += "\n# To enable the double, the minimum arch is 6.0\n";
                sContent += string.Format("target_compile_options({0} PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:-gencode arch=${1},code=${2}>)\n\n", m_sName, "{CUDA_CMP}", "{CUDA_SM}");

                return sContent;
            }
        }

        static string MakeCMake(EArch eArch)
        {
            string sContent = "cmake_minimum_required(VERSION 3.8 FATAL_ERROR)\n\n";

            sContent += "if (DEFINED NVCCROOT)\n";
            sContent += "    set(CMAKE_CUDA_COMPILER ${NVCCROOT})\n";
            sContent += "    MESSAGE(\"CMAKE_CUDA_COMPILER = ${CMAKE_CUDA_COMPILER}\")\n";
            sContent += "endif()\n\n";

            sContent += string.Format("set(CUDA_CMP \"{0}\")\n", ArchNames[(int)eArch]);
            sContent += string.Format("set(CUDA_SM \"{0}\")\n", CodeNames[(int)eArch]);

            sContent += "if (DEFINED CUDASM)\n";
            sContent += "    set(CUDA_CMP \"compute_${CUDASM}\")\n";
            sContent += "    set(CUDA_SM \"sm_${CUDASM}\")\n";
            sContent += "    if (NOT DEFINED CMAKE_CUDA_ARCHITECTURES)\n";
            sContent += "        set (CMAKE_CUDA_ARCHITECTURES ${CUDASM})\n";
            sContent += "    endif(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)\n";
            sContent += "endif()\n\n";
            sContent += "MESSAGE(\"Note: arch is ${CUDA_CMP} and ${CUDA_SM}.\")\n";
            sContent += "MESSAGE(\"52 for 970, 61 for GTX10, 70 for V100, 75 for RTX20, RTX16, 86 for RTX30\")\n";

            sContent += "project(QuantumLabProj LANGUAGES C CXX CUDA)\n\n";
            sContent += "set(CMAKE_GENERATOR_PLATFORM x64)\n\n";

            sContent += "# We start from CMAKE_SOURCE_DIR which should be /Code/CMake\n";
            sContent += "set(CMAKE_CURRENT_BINARY_DIR ${CMAKE_SOURCE_DIR}/../../Bin/Ubuntu)\n";
            sContent += "set(EXECUTABLE_OUTPUT_PATH  ${CMAKE_CURRENT_BINARY_DIR})\n";
            sContent += "set(LIBRARY_OUTPUT_PATH  ${CMAKE_CURRENT_BINARY_DIR})\n";

            sContent += "# This is our code file dir\n";
            sContent += "set(PROJECT_SOURCE_DIR ${CMAKE_SOURCE_DIR}/..)\n";

            sContent += "# Flags\n";
            sContent += "set(CMAKE_CUDA_FLAGS \"${CMAKE_CUDA_FLAGS} -O3\")\n";
            sContent += "set(CMAKE_CXX_FLAGS \"${CMAKE_CXX_FLAGS} -Ofast -Wall -Wno-unknown-pragmas -Wno-strict-overflow -Wno-class-memaccess\")\n";
            sContent += "add_definitions(-D_UBUNTU)\n";

            sContent += "MESSAGE(\"CMAKE_CUDA_FLAGS flag = ${CMAKE_CUDA_FLAGS}\")\n";
            sContent += "MESSAGE(\"CMAKE_CXX_FLAGS flag = ${CMAKE_CXX_FLAGS}\")\n\n";

            return sContent;
        }

        static void Main(string[] args)
        {
            //======================================
            //Config here
            string[][] sLibLinkWith = { new string[]{ }, new []{ "QuEST" } };
            string[] sLibNames = { "QuEST", "QuantumLab" };

            string[][] sAppLinkWith = { new[] { "QuantumLab" }, new[] { "QuantumLab" }, new[] { "QuantumLab" }, new[] { "QuantumLab" } };
            string[] sApplicationName = { "SimpleTest", "QKMeans", "SwapTest", "FermionSimulation" };
            EArch eDefaultArch = EArch.EArcSM61;

            string sContent = MakeCMake(eDefaultArch);
            //======================================
            for (int i = 0; i < sLibNames.Length; ++i)
            {
                CProjFile projApp = new CProjFile(sLibNames[i]);
                sContent += projApp.GetCMakeContentAsCudaLib(sLibLinkWith[i]);
            }

            for (int i = 0; i < sApplicationName.Length; ++i)
            {
                CProjFile projApp = new CProjFile(sApplicationName[i]);
                sContent += projApp.GetCMakeContentAsApplication(sAppLinkWith[i]);
            }

            sContent = sContent.Replace("\r\n", "\n");
            sContent = sContent.Replace("\n\r", "\n");
            sContent = sContent.Replace("\r", "\n");
            sContent = sContent.Replace("\\", "/");

            string sSolDir = Path.Combine(new[] { System.AppDomain.CurrentDomain.BaseDirectory, "../../" });
            File.WriteAllText(sSolDir + "Code/CMake/CMakeLists.txt", sContent);

            Console.WriteLine("work done, press enter to exit...");
            string byebye = Console.ReadLine();
        }
    }
}


