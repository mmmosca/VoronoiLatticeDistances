(*Copyright (C) 2018-2022,2023,2026 Marco M. Mosca

This file is part of VoronoiLatticeDistances.

VoronoiLatticeDistances is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

VoronoiLatticeDistances is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with VoronoiLatticeDistances. If not, see <https://www.gnu.org/licenses/>.
*)
open System.IO
open System.IO.Compression
open System.Net.Http
open System.Diagnostics
open System

let download (output_dir:string) (http_client: HttpClient) (uri: string) = 
    let path = Path.Combine(output_dir, Path.GetFileName(uri))
    async {
        let fileBytes = http_client.GetByteArrayAsync(uri).Result
        File.WriteAllBytes(path, fileBytes);
    } |> Async.RunSynchronously
    path

let start_process (proc: string) (args: string array) =
    let startInfo = ProcessStartInfo()
    for a in args do
        startInfo.ArgumentList.Add(a)
    startInfo.FileName <- proc
    startInfo.RedirectStandardOutput <- true
    startInfo.RedirectStandardError <- true
    startInfo.UseShellExecute <- false
    startInfo.CreateNoWindow <- true
    let proc = new Process(StartInfo = startInfo)
    proc.Start() |> ignore
    let stdout = proc.StandardOutput.ReadToEnd()
    let stderr = proc.StandardError.ReadToEnd()
    proc.WaitForExit()
    if proc.ExitCode <> 0 then
        failwithf "Process failed (%d): %s" proc.ExitCode stderr
    printfn "%s" stdout

let configure (src_dir:string) (variables: string array) = 
    let build_dir = Path.Combine(src_dir, "build")
    variables
    |> Array.append [|
        "-B"; build_dir; 
        "-S"; src_dir;
        "-L";
        sprintf @"-DCMAKE_INSTALL_PREFIX=%s" (Path.Combine(src_dir, "install-dir"));
        @"-DCMAKE_C_COMPILER=cl";
    |]
    |> start_process "cmake"
     
    build_dir

let install (build_type: string) (build_dir: string) =
    [|
        "--build"; build_dir;
        "--target"; "install";
        "--config"; build_type
    |]
    |> start_process "cmake"

let ctest (build_type: string) (build_dir: string) =
    [|
        "--test-dir"; build_dir;
        "--build-target"; "install";
        "--build-config"; build_type;
        "--extra-verbose"
    |]
    |> start_process "ctest"

let unzip (output_folder: string) (file: string) =
    ZipFile.ExtractToDirectory(
        file,
        output_folder,
        true
    )

let addToUserEnvPath (pathToFolder:string) =
    assert Directory.Exists(pathToFolder)
    let mutable path = Environment.GetEnvironmentVariable("Path", EnvironmentVariableTarget.User);
    if not (path.Contains(pathToFolder)) then
        path <- path + sprintf ";%s" pathToFolder;
        Environment.SetEnvironmentVariable("Path", path, EnvironmentVariableTarget.User);

let installDependencies ext_dir build_type =
    let http = new HttpClient()

    let src_nobuild_deps = [|
        (@"https://gitlab.com/libeigen/eigen/-/archive/3.3.7/eigen-3.3.7.zip", @"eigen-3.3.7");
        (@"https://github.com/project-gemmi/gemmi/archive/refs/tags/v0.3.3.zip", @"gemmi-0.3.3");
        (@"https://storage.googleapis.com/google-code-archive-downloads/v2/code.google.com/npackd/org.xmlsoft.LibXML-2.7.8.zip", @"libxml2-2.7.8.win32");
        (@"https://github.com/CGAL/cgal/releases/download/releases%2FCGAL-4.14.3/CGAL-4.14.3.zip", @"CGAL-4.14.3");
        (@"https://github.com/CGAL/cgal/releases/download/v5.2/CGAL-5.2-win64-auxiliary-libraries-gmp-mpfr.zip", @"auxiliary");
    |]
    let cgal_dir = Path.Combine(ext_dir, @"CGAL-4.14.3")
    for dep, folder_name in src_nobuild_deps do
        match Directory.Exists(Path.Combine(ext_dir,folder_name)) with
        | false ->
            printfn "%s - %s" folder_name "Downloading..."
            let path_to_file = download ext_dir http dep
            printfn "%s - %s" folder_name "Unzipping..."
            unzip ext_dir path_to_file
            File.Delete(path_to_file)
            match folder_name with
            | "gemmi-0.3.3" ->
                let small_hpp = Path.Combine(ext_dir, "gemmi-0.3.3/include/gemmi/small.hpp")
                assert File.Exists(small_hpp)
                let new_lines = [| 
                    @"  if (len > 1) {";
                    @"    if (label[len + 1] == '-')";
                    @"      dest->charge = -dest->charge;";
                    @"  }"
                |]
                let small_lines = 
                    System.IO.File.ReadAllLines(small_hpp)
                    |> Array.indexed 
                    |> Array.filter (fun (i,x) -> i <> 69 && i <> 70)
                    |> Array.map (fun (i,x) -> x) 
                    |> Array.insertManyAt 69 new_lines
                File.Delete(small_hpp)
                File.WriteAllLines(small_hpp, small_lines)
            | "libxml2-2.7.8.win32" ->
                Path.Combine(ext_dir, folder_name, @"bin") |> addToUserEnvPath 
            | "CGAL-4.14.3" ->
                let function_object_h = Path.Combine(cgal_dir, "include/CGAL/Cartesian/function_objects.h")
                assert File.Exists(function_object_h)
                let function_object_lines = 
                    System.IO.File.ReadAllLines(function_object_h)
                    |> Array.indexed 
                    |> Array.filter (fun (i,x) -> i <> 2209)
                    |> Array.map (fun (i,x) -> x) 
                File.Delete(function_object_h)
                File.WriteAllLines(function_object_h, function_object_lines)
                Environment.SetEnvironmentVariable("CGAL_DIR", cgal_dir, EnvironmentVariableTarget.User);
            | "auxiliary" ->
                Directory.Delete(Path.Combine(cgal_dir, folder_name, @"gmp"), true)
                Directory.Move(Path.Combine(ext_dir, folder_name, @"gmp"), Path.Combine(cgal_dir, folder_name, @"gmp"))
                Path.Combine(cgal_dir, @"auxiliary\gmp\lib") |> addToUserEnvPath
            | _ -> ()
        | true ->
            match folder_name with
            | "libxml2-2.7.8.win32" ->
                Path.Combine(ext_dir, folder_name, @"bin") |> addToUserEnvPath
            | "CGAL-4.14.3" ->
                Environment.SetEnvironmentVariable("CGAL_DIR", cgal_dir, EnvironmentVariableTarget.User);
            | "auxiliary" ->
                Path.Combine(cgal_dir, @"auxiliary\gmp\lib") |> addToUserEnvPath
            | _ -> ()


    let src_deps = [|
        (@"https://github.com/boostorg/boost/releases/download/boost-1.81.0/boost-1.81.0.zip", @"boost-1.81.0", [|@"-DBUILD_SHARED_LIBS=ON";@"-DBUILD_TESTING=OFF"|]);
        (@"https://www.zlib.net/zlib132.zip", @"zlib-1.3.2", [||]);
        (@"https://www.vtk.org/files/release/8.2/VTK-8.2.0.zip", @"VTK-8.2.0", [|
            @"-DBUILD_EXAMPLES=OFF";
            @"-DBUILD_TESTING=OFF";
            @"-DVTK_BUILD_ALL_MODULES=OFF";
            @"-DModule_vtkCommonCore=ON";
            @"-DModule_vtkCommonDataModel=ON";
            @"-DModule_vtkCommonColor=ON";
            @"-DModule_vtkCommonTransforms=ON";
            @"-DModule_vtkFiltersGeneral=ON";
            @"-DModule_vtkFiltersSources=ON";
            @"-DModule_vtkIOXML=ON";
            @"-DModule_vtkInteractionStyle=ON";
            @"-DModule_vtkRenderingCore=ON";
            @"-DModule_vtkRenderingOpenGL2=ON";
        |]);
    |]
    for dep, folder_name, variables in src_deps do
        match Directory.Exists(Path.Combine(ext_dir,folder_name)) with
        | false ->
            printfn "%s - %s" folder_name "Downloading..."
            let path_to_file = download ext_dir http dep
            printfn "%s - %s" folder_name "Unzipping..."
            unzip ext_dir path_to_file
            File.Delete(path_to_file)
            printfn "%s - %s" folder_name "Installing..."
            configure (Path.Combine(ext_dir, folder_name)) variables
            |> install build_type
            Path.Combine(ext_dir, folder_name, @"install-dir\bin")
            |> addToUserEnvPath

            match folder_name with
            | "boost-1.81.0" ->
                // Rename Boost library files for a naming issue in Boost 
                let bintorename = Path.Combine(ext_dir, folder_name, @"install-dir\bin") |> Directory.GetFiles
                let libtorename = Path.Combine(ext_dir, folder_name, @"install-dir\lib") |> Directory.GetFiles
                Array.append libtorename bintorename
                |> Array.iter (
                    fun f -> File.Copy(f, f.Replace(@"vc144",@"vc143"))
                )
            | _ -> ()
        | true ->
            Path.Combine(ext_dir, folder_name, @"install-dir\bin")
            |> addToUserEnvPath

[<EntryPoint>]
let main args =
    let subcommand = args[0]
    let ext_dir = Path.GetFullPath(@"..\External")
    let build_dir = @"..\build"

    if Directory.Exists build_dir then
        Directory.Delete(build_dir, true)
    match subcommand with
    | "installdeps" ->
        let build_type = args[1]
        match build_type with
        | "Debug" | "Release" -> 
            installDependencies ext_dir build_type
        | _ ->
            failwithf "Unknown build type: %s" subcommand
    | "install" ->
        let build_type = args[1]
        match build_type with
        | "Debug" ->
            installDependencies ext_dir build_type
            printfn "%s - %s" "Project Voronoi" "Installing..."
            configure (Path.GetFullPath(@"..\")) [|@"-DSHARED=OFF";@"-DBUILD_TESTING=ON"|] |> install build_type
        | "Release" -> 
            installDependencies ext_dir build_type
            printfn "%s - %s" "Project Voronoi" "Installing..."
            configure (Path.GetFullPath(@"..\")) [|@"-DSHARED=ON";@"-DBUILD_TESTING=OFF"|] |> install build_type
        | _ ->
            failwithf "Unknown build type: %s" subcommand
        Path.GetFullPath(@"..\install-dir\bin") |> addToUserEnvPath
        printfn "%s - %s" "Project Voronoi" "Completed!"
        printfn "%s - %s" "Project Voronoi" "Restart your command prompt to use the executables!"
    | "tests" ->
        let build_type = "Debug"
        installDependencies ext_dir build_type
        configure (Path.GetFullPath(@"..\")) [|@"-DSHARED=OFF";@"-DBUILD_TESTING=ON"|]
        |> fun bdir ->
            install build_type bdir
            ctest build_type bdir
    | _ -> 
        failwithf "Unknown command: %s" subcommand

    0
