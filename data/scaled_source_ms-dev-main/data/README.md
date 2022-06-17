# DATAPREP_sattag_processing

Basic sattag processing for the most commonly requested streams and products.

- status
- locations
- behavior
- series / series range

main products are merging in goniometer messages and some basic filtering including our rather stringent pressure transducer health criteria.

**note:** you'll notice that fastloc isn't available at the moment here. All the portal fastloc is available in tag\_datastreams\_raw, but it has not been combined with goniometer fastloc due to technical difficulties. If fastloc data is a priority for you please make an issue and send me an issue and I'll move it up the list.

see below for details

## HOW TO USE
**WARNING:**
I've done my best to QAQC, but please double check to make sure things make sense. I've done only the most basic error detection and filtering here because I do not know what your downstream analysis is. Please institute your own filtering and tests. If you find what you think is an error or have any doubts please file and issue and email wrc14 [at] duke [dot] edu.

This repo is not meant to be the primary method for sharing the data, but rather as a way to trouble shoot and develop the code pipeline for processing.

for accessing the data you probably want to look at:
- tag\_datastreams\_processed: https://duke.box.com/s/26yxi6zxnoz801vt0vze4x9ynj0b1rxd
- tag\_datastreams\_raw: https://duke.box.com/s/0xq7irgf8w0h3dr10rn9gv8gp5avia49
- tagdeploy\_metadata: https://duke.box.com/s/nd0ypwzvpnp7kp0o1ipqdp8fuunkmm4k
- cee\_metadata\_flat: https://duke.box.com/s/mi9ravh1nvkq04xjcq4pnx476epeulj1
- ATLANTIC-BRS_CEE Metadata: https://duke.box.com/s/vtiawx7xznhporib8s1zjlef3ymu14on
- goniometer hit locations / strengths:  https://duke.box.com/s/puxd3shhmaz6hffywynv5n6xxgf2s4hu 


## DETAILS:


## STATUS
processed status message data for all BRS sattags

Notes:
- If no status messages ever arrived then no file will be present here.
- If all status messages are filtered out, then a file with no rows will be present here.

### L0\_raw
Not strictly speaking raw and differs in several important ways from for instance the data streams downloaded from the wildlife computers portal (or from the raw tag archives). Also lightly filtered. See below.

Differences to wildlife computers portal data:
- Dates have been converted into datenum (seconds) since the UNIX epoch 1970-01-01. Excel is horrible with dates (especially in CSVs) so this is a more safe / portable format. If it is annoying I can put strings back in.
- DeployIDs have been corrected where necessary. There is no naming conflict with DeployIDs within the BRS, but there is between projects. The presence or absence of the '\_DUML' suffix is meaningful.
- columns have been added for forward compatibility where necessary. These include `Tilt` (added first), `S11`, `ReleaseWetDry`, and `ReleaseTemperature` (added later). These are all always NA for our tags so have been populated with those values on older tags.

Filtering:
- Any status messages before the deploy date have been filtered out. Sometimes tags are tested and the messages don't always get removed.
- Technically a filter was applied to cut off any spurious messages after our tag had most likely stopped transmitting although I don't think this actually applies to the status stream currently.

### L1_crc
Filtered only for messages that pass the cyclic redundancy check. Occasionally messages fail the CRC, but more often for some reason no CRC is able to be performed. I generally take a conservative approach and only accept CRC'd messages for identifying pressure transducer failures, for example.

Otherwise these tables should be identical to L0.

### L2_crc_gonio
Same as L1, but merged in with any status messages received on the goniometer.

Additional columns:

- `original` indicates whether a given line is from either 'portal' (i.e., the wildlife computers data portal) or 'gonio' (i.e., received on the goniometer)

Notes:

- duplicates have been removed, but it is somewhat random whether the gonio or portal line is removed, so the `original` column is not useful for assessing satellite performance for instance. This can be achieved however by comparing the separate 'gonio' and 'portal' data directories found in the raw tag archive.

## LOCATIONS
processed locations message data for all BRS sattags.

Notes:

- If no location messages ever arrived then no file will be present here.
- fastloc positions are *NOT* included (see also issue #2).

### L0_raw
Not strictly speaking raw and differs in several important ways from for instance the data streams downloaded from the wildlife computers portal (or from the raw tag archives). Also lightly filtered. See below.

Differences to wildlife computers portal data:

- Dates have been converted into datenum (seconds) since the UNIX epoch 1970-01-01. Excel is horrible with dates (especially in CSVs) so this is a more safe / portable format. If it is annoying I can put strings back in.
- DeployIDs have been corrected where neccessary. There is no naming conflict with DeployIDs within the BRS, but there is between projects. The presence or abence of the '\_DUML' suffix is meaningful.
- an `originaldate` column has been added with the format `%H:%M:%S %d-%b-%Y` (in the R POSIX functions format convention). This is a funny format, but is often used in the wildlife computers spreadsheets

Filtering:

- Any location message before the deploy date have been filtered out. Sometimes tags are tested and the messages don't always get removed.
- Technically a filter was applied to cut off any spurious messages after our tag had most likely stopped transmitting although I don't think this actually applies to the locations stream currently.

### L1_kalman (PROVISIONAL SEE NOTE BELOW)
Tags prior to 2019 had their locations processed using the old argos least squares solver. These tags were reprocessed by cls for us to produce new positions using the kalman filter algorithm. Tags starting in 2019 were processed using the kalman filter from the beginning and for these tags this L1 is identical to L2.

For tags prior to 2019, the cls kalman position output was simply reformatted to appear similar to the wildlife computers portal spreadsheet output. There are some small differences worth noting. I'll simply call the wildlife computers output 'locations' and the cls output 'kalman', but note that they are both kalman filtered in this case.

- locations has a field `Error.Ellipse.orientation`, but kalman does not
- kalman includes fields that locations does not (and therefore are not included in L1). Although these are still available in the raw tag archive (in the kalman directory) if they are desired. Tags in 2019 and later do not have a kalman directory in their raw tag archive, but many of the field are available in the `argos`, `rawargos`, or `all` data streams.
  - satellite id 
  - messages 
  - pass 
  - frequency
  - GDOP
  - Dist. Subsat Track
  - second set (solve) of longitude and latitude (tend to be identical in the kalman datastream)

<b>IMPORTANT NOTE:</b>

There is no easy 1:1 correspondence between the `argos` datastream for instance from the wildlife computers data portal spreadsheet and the kalman filtered spreadsheets we received from cls. While times, satellites, qualities, all can be very similar for many rows, there are a substantial number of rows that do not match between these two data sources. In addition, there are small time discrepancies (on the order of seconds) between rows that do appear to match. cls has been unable to provide a detailed explanation of the differences. wildlife computers notes that any updates made on cls' side should propagate through the system within 48 hours and be reflected in the wildlife computers dataportal spreadsheet outputs. 

For now, exercise caution before using L1 and do your own double checking!

## BEHAVIOR
processed behavior message data for all BRS sattags

Notes:

- If no behavior messages ever arrived then no file will be present here.
- If all behavior messages are filtered out, then a file with no rows will be present here.

### L0_raw
Not strictly speaking raw and differs in several important ways from for instance the data streams downloaded from the wildlife computers portal (or from the raw tag archives). Also lightly filtered. See below.

Differences to wildlife computers portal data:

- Dates have been converted into datenum (seconds) since the UNIX epoch 1970-01-01. Excel is horrible with dates (especially in CSVs) so this is a more safe / portable format. If it is annoying I can put strings back in.
- DeployIDs have been corrected where necessary. There is no naming conflict with DeployIDs within the BRS, but there is between projects. The presence or absence of the '\_DUML' suffix is meaningful.

Filtering:

- Any behavior message before the deploy date have been filtered out. Sometimes tags are tested and the messages don't always get removed.
- A cut off was applied to remove spurious messages after our tag had most likely stopped transmitting

### L1_gonio
Same as L0, but merged with any behavior messages received on the goniometer.

Additional columns:

- `original` indicates whether a given line is from either 'portal' (i.e., the wildlife computers data portal) or 'gonio' (i.e., received on the goniometer)

Notes:

- duplicates have been removed, but it is somewhat random whether the gonio or portal line is removed, so the `original` column is not useful for assessing satellite performance for instance. This can be achieved however by comparing the separate 'gonio' and 'portal' data directories found in the raw tag archive.

### L2_gonio_pressuresensor
Same as L1, but with an additional filter to try to identify the most obvious pressure transducer failures. 

A record is flagged if 2 status messages (CRC'd, combied with goniometer received messages) were logged where the |zerodepth| > 10 (meters). The cutoff time is rolled back to the last good status message (|zerodepth| <= 10) prior to the 1st of the bad messages.

Otherwise the cutoff is the latest CRC'd (good) status message.

This is actually a rather severe filtering step, especially for tags which did not successfully send many status messages.

As noted above, if all messages were filtered out there will be an empty file.

Also Note: when there is truncation, the End time of the last message (beh$What == "Message") is not changed. I did this so you could see what it was supposed to be, but if it is a problem it is easy to update.

## SERIES / SERIESRANGE
processed series and seriesrange message data for all BRS sattags

Notes:

- If no behavior messages ever arrived then no file will be present here.
- If all behavior messages are filtered out, then a file with no rows will be present here.

### L0_raw

Not strictly speaking raw and differs in several important ways from for instance the data streams downloaded from the wildlife computers portal (or from the raw tag archives). Also lightly filtered. See below.

Differences to wildlife computers portal data:

- Dates have been converted into datenum (seconds) since the UNIX epoch 1970-01-01. Excel is horrible with dates (especially in CSVs) so this is a more safe / portable format. If it is annoying I can put strings back in.
  - datenums are in a new column called `Date` in series and in the existing `Start` and `End` in seriesrange. In series the original dates are split between day and time and those are left in their original formats
- DeployIDs have been corrected where necessary. There is no naming conflict with DeployIDs within the BRS, but there is between projects. The presence or absence of the '\_DUML' suffix is meaningful.

Filtering:

- Any behavior message before the deploy date have been filtered out. Sometimes tags are tested and the messages don't always get removed.
- A cut off was applied to remove spurious messages after our tag had most likely stopped transmitting

### L1_gonio

Same as L0, but merged with any series and seriesrange messages received on the goniometer.

Additional columns:

- `original` indicates whether a given line is from either 'portal' (i.e., the wildlife computers data portal) or 'gonio' (i.e., received on the goniometer)

Notes:

- duplicates have been removed, but it is somewhat random whether the gonio or portal line is removed, so the `original` column is not useful for assessing satellite performance for instance. This can be achieved however by comparing the separate 'gonio' and 'portal' data directories found in the raw tag archive.

### L2_gonio_pressuresensor

Same as L1, but with an additional filter to try to identify the most obvious pressure transducer failures. 

rules:
1. flagged if 2 CRC'd status messages with |zerodepth| > 10
2. if all drift is in the shallow direction cutoff to the last good sta
3. if all drift is in the deeper direction, truncate to where the 8 hour rolling min in the series falls below 20 meters (calculated from all good baseline)
4. if there is a mix of deeper and shallow drift, truncate to the last negative drift with less than 20 m rolling min, or the last good status message, whichever is later.

As a quick justification for this, it is possible to detect negative (deeper) drift by looking at the minimum bin depths are falling into, but not positive drift (shallower) since there is no bin above 0 m due to constraints of the tag's depth encoding. This is still quite conservative, but since the series data does give a hint at negative drift we use that to identify when the onset of drift occured between the last good status message and the first bad (where possible).

As noted above, if all messages were filtered out there will be an empty file.

Also Note: when there is truncation, the End time of the last SeriesRange message is not changed. I did this so you could see what it was supposed to be, but if it is a problem it is easy to update.
