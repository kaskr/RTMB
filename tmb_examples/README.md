## RTMB examples

Examples in this folder are translations of same examples from the [TMB examples folder](https://github.com/kaskr/adcomp/tree/master/tmb_examples).

## Correctness checks

Assuming that projects RTMB and adcomp are cloned side by side.
Start by copying expected output from adcomp:

```shell
make expected_output
```

Then run tests using:

```shell
make
make report
cat REPORT.md
```
