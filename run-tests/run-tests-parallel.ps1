param(
    [Parameter(Mandatory=$true)][string]$deborah,
    [switch]$saveOutputs = $false
)

$scenarios = Get-ChildItem -Recurse -Path "../scenarios/" -Filter "*.conf";

$commands = @(
    "lb",
    "lb --lb-experimental",
    "ludb",
    "ludb --ludb-nnested-sta",
    "foi-output-curve"
);

$outfilesFolder = "last-test-outputs";
if($saveOutputs)
{
    if(Test-Path($outfilesFolder))
    {
        Get-ChildItem -Recurse -Path $outfilesFolder | Remove-Item -Recurse;
    }
    else
    {
        New-Item -ItemType Directory -Path $outfilesFolder;
    }

    foreach($command in $commands)
    {
        New-Item -ItemType Directory -Path "$outfilesFolder/$command";
    }
}

function getOutfile($command, [System.IO.FileInfo] $scenario)
{
    if($saveOutputs)
    {
        return "$outfilesFolder/$command/$($scenario.BaseName).log"
    }
    else
    {
        return New-TemporaryFile;
    }    
}

$tests = @();
$currentId = 0;

foreach($command in $commands)
{
    foreach($scenario in $scenarios)
    {
        $tests += @{
            Id = $currentId
            deborah = $deborah
            command = $command
            scenario = $scenario
            outfile = getOutfile $command $scenario
        }
        $currentId += 1;
    }
}

# Create a hashtable for process.
# Keys should be ID's of the processes
$origin = @{}
$tests | Foreach-Object {$origin.($_.Id) = @{}}

# Create synced hashtable
$sync = [System.Collections.Hashtable]::Synchronized($origin)


$stopwatch = [System.Diagnostics.Stopwatch]::StartNew();

$job = $tests | ForEach-Object -AsJob -Parallel {
    $syncCopy = $using:sync;
    $process = $syncCopy.$($_.Id);

    # Initialize process report
    $process.Id = $_.Id
    $process.Completed = $false;
    $process.Log = @();
    $process.Success = $false;
    $process.Time = "";

    # Run test
    function AreFilesEqual($a, $b)
    {
        $hashA = (Get-FileHash $a).Hash;
        $hashB = (Get-FileHash $b).Hash;

        if($hashA -eq $hashB) {
            return $true;
        }
        else {
            return $false;
        }
    }

    $deborah = $_.deborah;
    $command = $_.command;
    $scenario = $_.scenario;
    $outfile = $_.outfile;

    $process.Title = "$command $scenario";

    $process.Log += "Testing $command on $($scenario.Name)";

    $process.Log += "Run output: $outfile";
    $refFile = "$command/$($scenario.BaseName).log";
    $process.Log += "Reference file: $refFile";

    $expression = "$deborah $scenario --det-output --$command > `"$outfile`"";
    $process.Log += "Command: $expression";

    $stopwatch = [System.Diagnostics.Stopwatch]::StartNew();
    Invoke-Expression $expression;
    $stopwatch.Stop();

    $process.Time = $stopwatch.Elapsed;

    $passed = AreFilesEqual $outfile $refFile;
    if(-Not $passed)
    {
        $process.Log += "TEST FAILED";
        $process.Log += "Reference is left, test output is right";
        $process.Log += Compare-Object -ReferenceObject (Get-Content $refFile) -DifferenceObject (Get-Content $outfile);

        $process.Success = $false;
    }
    else {
        $process.Log += "TEST PASSED";
        $process.Success = $true;
    }

    # Mark process as completed
    $process.Completed = $true;
};

while($job.State -eq 'Running')
{
    $completed = $sync.Keys |
        # Where-Object { -not [string]::IsNullOrEmpty($sync.$_.keys) } |
        Where-Object { $sync.$_.Completed };

    $runsLeft = $currentId - $completed.Count;
    $timeLeft = $stopwatch.Elapsed.ToString('dd\.hh\:mm\:ss');

    if($runsLeft -eq 1)
    {
        $title = $sync[($sync.Keys | Where-Object { -not $sync.$_.Completed })[0]].Title;
        $activity = "1 run left: $title $timeLeft";
    }
    else
    {
        $activity = "$runsLeft runs left. $timeLeft";
    }

    Write-Progress -Activity $activity -PercentComplete $($completed.Count * 100 / $currentId)

    # $sync.Keys | Foreach-Object {
    #     # If key is not defined, ignore
    #     if(![string]::IsNullOrEmpty($sync.$_.keys))
    #     {
    #         # Create parameter hashtable to splat
    #         $param = $sync.$_

    #         # Execute Write-Progress
    #         Write-Progress @param
    #     }
    # }

    # Wait to refresh to not overload gui
    Start-Sleep -Seconds 1
}

$stopwatch.Stop();

$results = @();
foreach($key in $origin.Keys)
{
    $results += $origin[$key];
}

$totalTests = $results.Length;
$passed = ( $results | Where-Object { $_.Success }).count;

Write-Output "Passed $passed out of $totalTests tests, in $($stopwatch.Elapsed)";

foreach($result in $results)
{
    if($result.Success)
    {
        if($Verbose)
        {
            Write-Output $result.Log
        }
    }
    else
    {
        Write-Output $result.Log
    }
}