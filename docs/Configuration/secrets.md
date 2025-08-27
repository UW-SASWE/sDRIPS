# Secrets Configuration File
Sensitive credentials, such as [***PPS credentials***](/en/latest/Quick_Start/Getting_ready/#pps-registration) or [***Google Earth Engine***](/en/latest/Quick_Start/Getting_ready/#setup-gee) service account information, are stored in a dedicated secrets configuration file (`secrets.yaml`). Separating this information from the main workflow enhances **security**, prevents accidental exposure of sensitive data, and improves **portability** across systems and users.

```yaml
GEE_Account:
  username: <GEE_USERNAME> # Google Earth Engine service account
  key_file: <GEE_KEY_FILE> # Path to the key file

IMERG_Account:
  username: <IMERG_USERNAME> # IMERG username
  password: <IMERG_PASSWORD> # IMERG password
```
!!! tip_note "Tip"
    Users should populate the secrets file immediately after completing the [***initialization***](/en/latest/Quick_Start/Installation/#initialization). Doing so allows sDRIPS to run its [***test suite***](Quick_Start/Installation/#testing), which verifies both the proper installation of the software and the correctness of the credentials stored in the secrets file.