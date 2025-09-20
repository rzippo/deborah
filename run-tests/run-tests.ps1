#!/bin/pwsh

param(
    [switch]$continueOnFail = $false,
    [switch]$verbose = $false,
    [string]$deborah = "deborah"
)

$scenarios = Get-ChildItem -Recurse -Path "../scenarios/" -Filter "*.conf";

$commands = @(
    "lb",
    "lb --lb-experimental",
    "ludb",
    "ludb --ludb-nnested-sta",
    "foi-output-curve"
);

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

$stopwatch = [System.Diagnostics.Stopwatch]::StartNew();

$nTests = 0;
$nTestsSuccess = 0;
$nTestsSkipped = 0;
$nTestsFailed = 0;

foreach($command in $commands)
{
    foreach($scenario in $scenarios)
    {
        $nTests++;
        Write-Output "Testing $command on $($scenario.Name)";

        $refFile = "$command/$($scenario.BaseName).log";
        if(-not(Test-Path -Path $refFile -PathType Leaf)) {
            Write-Output "TEST SKIPPED: ref file does not exist";
            $nTestsSkipped++;
            continue;
        }

        $tmpFile = New-TemporaryFile;
        if($verbose)
        {
            Write-Output "Run output: $tmpFile";
        }

        $expression = "$deborah $scenario --det-output --$command > $tmpFile";
        if($verbose)
        {
            Write-Output "Command: $expression";
        }
        Invoke-Expression $expression;

        $passed = AreFilesEqual $tmpFile $refFile;
        if(-Not $passed)
        {
            $nTestsFailed++;
            Write-Output "TEST FAILED";
            Write-Output "Reference is left, test output is right";
            Compare-Object -ReferenceObject (Get-Content $refFile) -DifferenceObject (Get-Content $tmpFile)

            if(-Not $continueOnFail)
            {
                exit;
            }
        }
        else {
            $nTestsSuccess++;
            Write-Output "TEST PASSED"
        }
    }    
}

$stopwatch.Stop();
Write-Output "Running time: $($stopwatch.Elapsed)";
Write-Output "$nTests tests: $nTestsSuccess passed, $nTestsSkipped skipped, $nTestsFailed failed";