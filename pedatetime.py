"""
pedatetime.py
-----------------

A lightweight, custom datetime and timedelta implementation intended for
astronomical and historical time calculations.

Key characteristics:
- No timezone support (all calculations are strictly naive).
- Uses second-by-second arithmetic to enssure deterministic historical behavior.
- Handles the Gregorian calendar cutover (skips 1582-10-05 → 1582-10-14).
- Provides custom normalization of fields (seconds, minutes, hours, days, months).
- Includes a simple timedelta class based on total seconds.

This module is designed for astronomical purposes where precise
control over historical date transitions is required, and should NOT be used
as a replacement for Python's built-in datetime library for modern applications,
especially those requiring timezone or leap second support.

"""


# ---------------------------------------------------------------------------
# Timedelta implementation
# ---------------------------------------------------------------------------
class timedelta:
    """
    Represents a duration expressed in days, hours, minutes, and seconds.
    Internally stores the value as total seconds (may be negative).
    """

    def __init__(self, days: int, hours: int, minutes: int, seconds: int) -> None:
        # Type validation
        if not all(isinstance(v, int) for v in (days, hours, minutes, seconds)):
            raise TypeError("timedelta arguments must be integers.")

        # Convert everything into total seconds
        self.total_seconds = days * 86400 + hours * 3600 + minutes * 60 + seconds

    def __add__(self, other: "timedelta") -> "timedelta":
        if not isinstance(other, timedelta):
            return NotImplemented
        return timedelta(0, 0, 0, self.total_seconds + other.total_seconds)

    def __sub__(self, other: "timedelta") -> "timedelta":
        if not isinstance(other, timedelta):
            return NotImplemented
        return timedelta(0, 0, 0, self.total_seconds - other.total_seconds)

    def __neg__(self) -> "timedelta":
        return timedelta(0, 0, 0, -self.total_seconds)

    # Comparison operators
    def __eq__(self, other: object) -> bool:
        if not isinstance(other, timedelta):
            return NotImplemented
        return self.total_seconds == other.total_seconds

    def __lt__(self, other: "timedelta") -> bool:
        return self.total_seconds < other.total_seconds

    def __le__(self, other: "timedelta") -> bool:
        return self.total_seconds <= other.total_seconds

    def __gt__(self, other: "timedelta") -> bool:
        return self.total_seconds > other.total_seconds

    def __ge__(self, other: "timedelta") -> bool:
        return self.total_seconds >= other.total_seconds


# ---------------------------------------------------------------------------
# Helper function: days in a given month
# ---------------------------------------------------------------------------
def max_day_in_month(year: int, month: int) -> int:
    """
    Returns the number of days in a given month of a given year.
    Handles leap years. Returns None if the month is invalid.
    """
    if not isinstance(year, int) or not isinstance(month, int):
        raise TypeError("year and month must be integers")

    if month < 1 or month > 12:
        return 0

    month_lengths = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

    # Determine leap year
    leap = year % 4 == 0 and (year % 100 != 0 or year % 400 == 0)
    if leap:
        month_lengths[1] = 29

    return month_lengths[month - 1]


# ---------------------------------------------------------------------------
# Datetime implementation
# ---------------------------------------------------------------------------
class datetime:
    """
    Custom datetime class supporting:
    - manual normalization
    - Gregorian cutover adjustment
    - second-by-second arithmetic
    """

    def __init__(
        self,
        yr: int,
        mn: int,
        dy: int,
        hh: int,
        mm: int,
        ss: int,
    ) -> None:
        # Type checking
        if not all(isinstance(v, int) for v in (yr, mn, dy, hh, mm, ss)):
            raise TypeError("All datetime components must be integers.")

        self.year = yr
        self.month = mn
        self.day = dy
        self.hour = hh
        self.minute = mm
        self.second = ss

    # ----------------------------------------------------------------------
    # Normalization
    # ----------------------------------------------------------------------
    def check(self) -> None:
        """
        Normalizes the datetime fields:
        - Ensures seconds, minutes, hours stay in valid ranges.
        - Adjusts days and months with overflow/underflow.
        """
        # Normalize seconds
        while self.second < 0:
            self.second += 60
            self.minute -= 1
        while self.second > 59:
            self.second -= 60
            self.minute += 1

        # Normalize minutes
        while self.minute < 0:
            self.minute += 60
            self.hour -= 1
        while self.minute > 59:
            self.minute -= 60
            self.hour += 1

        # Normalize hours
        while self.hour < 0:
            self.hour += 24
            self.day -= 1
        while self.hour > 23:
            self.hour -= 24
            self.day += 1

        # Normalize days backward
        while self.day < 1:
            self.month -= 1
            if self.month < 1:
                self.month = 12
                self.year -= 1
            self.day += max_day_in_month(self.year, self.month)

        # Normalize days forward
        while self.day > max_day_in_month(self.year, self.month):
            self.day -= max_day_in_month(self.year, self.month)
            self.month += 1
            if self.month > 12:
                self.month = 1
                self.year += 1

        # Normalize months
        while self.month < 1:
            self.month += 12
            self.year -= 1
        while self.month > 12:
            self.month -= 12
            self.year += 1

    # ----------------------------------------------------------------------
    # Gregorian correction
    # ----------------------------------------------------------------------
    def checkGreg(self, addsub: bool) -> None:
        """
        Adjusts dates around the Gregorian cutover (1582-10-05 → 1582-10-14).
        `addsub=True` means addition; False means subtraction.
        """
        if addsub:
            if self.year == 1582 and self.month == 10 and 5 <= self.day <= 14:
                self.day = 15

        else:
            if self.year == 1582 and self.month == 10 and 5 <= self.day <= 14:
                self.day = 4

    # ----------------------------------------------------------------------
    # Increment operations
    # ----------------------------------------------------------------------
    def add_second(self) -> None:
        self.second += 1
        self.check()
        self.checkGreg(True)

    def add_minute(self) -> None:
        self.minute += 1
        self.check()
        self.checkGreg(True)

    def add_hour(self) -> None:
        self.hour += 1
        self.check()
        self.checkGreg(True)

    def add_day(self) -> None:
        self.day += 1
        self.check()
        self.checkGreg(True)

    # ----------------------------------------------------------------------
    # Decrement operations
    # ----------------------------------------------------------------------
    def sub_second(self) -> None:
        self.second -= 1
        self.check()
        self.checkGreg(False)

    def sub_minute(self) -> None:
        self.minute -= 1
        self.check()
        self.checkGreg(False)

    def sub_hour(self) -> None:
        self.hour -= 1
        self.check()
        self.checkGreg(False)

    def sub_day(self) -> None:
        self.day -= 1
        self.check()
        self.checkGreg(True)

    # ----------------------------------------------------------------------
    # Formatting
    # ----------------------------------------------------------------------
    def isoformat(self) -> str:
        """Return an ISO-8601 formatted datetime string."""
        return (
            f"{self.year:04d}-{self.month:02d}-{self.day:02d}"
            f"T{self.hour:02d}:{self.minute:02d}:{self.second:02d}"
        )

    # ----------------------------------------------------------------------
    # Comparisons
    # ----------------------------------------------------------------------
    def __eq__(self, other: object) -> bool:
        if not isinstance(other, datetime):
            return NotImplemented
        return (
            self.year,
            self.month,
            self.day,
            self.hour,
            self.minute,
            self.second,
        ) == (
            other.year,
            other.month,
            other.day,
            other.hour,
            other.minute,
            other.second,
        )

    def __lt__(self, other: "datetime") -> bool:
        return (
            self.year,
            self.month,
            self.day,
            self.hour,
            self.minute,
            self.second,
        ) < (other.year, other.month, other.day, other.hour, other.minute, other.second)

    def __le__(self, other: "datetime") -> bool:
        return self == other or self < other

    def __gt__(self, other: "datetime") -> bool:
        return not self <= other

    def __ge__(self, other: "datetime") -> bool:
        return not self < other

    def __ne__(self, other: object) -> bool:
        return not self == other

    # ----------------------------------------------------------------------
    # Utilities
    # ----------------------------------------------------------------------
    def copy(self) -> "datetime":
        """Returns a fresh copy of the datetime instance."""
        return datetime(
            self.year, self.month, self.day, self.hour, self.minute, self.second
        )

    # Repeated incremental operations (slow but preserves original behavior)
    def add_seconds(self, s: int) -> None:
        for _ in range(s):
            self.add_second()

    def add_minutes(self, m: int) -> None:
        for _ in range(m):
            self.add_minute()

    def add_hours(self, h: int) -> None:
        for _ in range(h):
            self.add_hour()

    def add_days(self, d: int) -> None:
        for _ in range(d):
            self.add_day()

    def sub_seconds(self, s: int) -> None:
        for _ in range(s):
            self.sub_second()

    def sub_minutes(self, m: int) -> None:
        for _ in range(m):
            self.sub_minute()

    def sub_hours(self, h: int) -> None:
        for _ in range(h):
            self.sub_hour()

    def sub_days(self, d: int) -> None:
        for _ in range(d):
            self.sub_day()

    # ----------------------------------------------------------------------
    # Arithmetic
    # ----------------------------------------------------------------------
    def __add__(self, td: timedelta) -> "datetime":
        if not isinstance(td, timedelta):
            return NotImplemented

        result = self.copy()

        # Apply the timedelta one second at a time
        if td.total_seconds >= 0:
            result.add_seconds(td.total_seconds)
        else:
            result.sub_seconds(-td.total_seconds)

        return result

    def __sub__(self, other: object):
        # datetime - timedelta
        if isinstance(other, timedelta):
            return self + (-other)

        # datetime - datetime => timedelta
        if isinstance(other, datetime):
            a = other.copy()
            b = self.copy()
            seconds = 0

            # Step second-by-second to compute delta
            if b >= a:
                while a < b:
                    a.add_second()
                    seconds += 1
            else:
                while a > b:
                    a.sub_second()
                    seconds -= 1

            return timedelta(0, 0, 0, seconds)

        return NotImplemented
