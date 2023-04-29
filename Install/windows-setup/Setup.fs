(*
 * Permission is granted to copy, distribute and/or modify the documents
 * in this directory and its subdirectories unless otherwise stated under
 * the terms of the GNU Free Documentation License, Version 1.1 or any later version 
 * published by the Free Software Foundation; with no Invariant Sections, 
 * no Front-Cover Texts and no Back-Cover Texts. A copy of the license 
 * is available at the website of the GNU Project.
 * The programs and code snippets in this directory and its subdirectories
 * are free software; you can redistribute them and/or modify it under the 
 * terms of the GNU General Public License as published by the Free Software 
 * Foundation; either version 2 of the License, or (at your option) any later
 * version. This code is distributed in the hope that it will be useful, 
 * but WITHOUT ANY WARRANTY; without even the implied warranty of 
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 * 
 * Author Marco M. Mosca, email: marcomichele.mosca@gmail.com
*)
open System.IO
open System.Net.Http
open System.Diagnostics

let unzip (output_folder: string) (file: string) =
    let startInfo = ProcessStartInfo()
    startInfo.FileName <- "7z.exe"
    for a in [| "x"; "-tzip"; file; sprintf "-o%s" output_folder |] do
        startInfo.ArgumentList.Add(a)
    startInfo.RedirectStandardOutput <- false
    startInfo.RedirectStandardError <- false
    startInfo.UseShellExecute <- false
    startInfo.CreateNoWindow <- true
    use p = new Process()
    p.StartInfo <- startInfo
    p.Start() |> ignore
    p.WaitForExit()

[<EntryPoint>]
let main args =
    
    let ext_dir = "../../External"
    assert Directory.Exists(ext_dir)
    let src_deps = [|
        @"https://gitlab.com/libeigen/eigen/-/archive/3.3.7/eigen-3.3.7.zip";
        @"https://github.com/CGAL/cgal/releases/download/releases%2FCGAL-4.14.3/CGAL-4.14.3.zip";
        @"https://github.com/project-gemmi/gemmi/archive/refs/tags/v0.3.3.zip";
    |]

    let http = new HttpClient()
    for dep in src_deps do
        async {
            let path = Path.Combine(ext_dir, Path.GetFileName(dep))
            printfn "Downloading and unzipping %s..." dep
            let fileBytes = http.GetByteArrayAsync(dep).Result
            File.WriteAllBytes(path, fileBytes);
            unzip ext_dir path
            File.Delete(path)
        } |> Async.RunSynchronously

    Directory.Move(Path.Combine(ext_dir,"eigen-3.3.7"), Path.Combine(ext_dir,"Eigen-3.3.7"))
    Directory.Move(Path.Combine(ext_dir,"gemmi-0.3.3"), Path.Combine(ext_dir,"Gemmi"))

    let function_object_h = Path.Combine(ext_dir, "CGAL-4.14.3/include/CGAL/Cartesian/function_objects.h")
    assert File.Exists(function_object_h)
    let function_object_lines = 
        System.IO.File.ReadAllLines(function_object_h)
        |> Array.indexed 
        |> Array.filter (fun (i,x) -> i <> 2209)
        |> Array.map (fun (i,x) -> x) 
    File.Delete(function_object_h)
    File.WriteAllLines(function_object_h, function_object_lines)

    let small_hpp = Path.Combine(ext_dir, "Gemmi/include/gemmi/small.hpp")
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

    0